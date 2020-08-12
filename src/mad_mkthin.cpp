/* mad_mkthin.cpp

 Thick to thin lens converter makethin. Helmut Burkhardt

 Major steps
 2001, 2002 early versions by Mark Hayes
 2005 : Standard selection SELECT,FLAG=makethin,RANGE=range,CLASS=class,PATTERN=pattern[,FULL][,CLEAR]; Implementation of slicing for solenoids
 2012 : Extension of TEAPOT slicing to n>4
 2013 : Keep thick elements if slice number <1, code now C++, thick slicing for quadrupoles, automatic generation of dipedge elements for dipoles
 2014 : Thick bend slicing, with or without dipedge
 2016 : Extra markers, sbend_from_rbend transmit aper_tol
 2018 : Thick solenoid slicing, write bend angle to multipole if different from k0*l
 2019 : New elements from element definition, all attributes enabled
 2020 : Tapering, wire compensation

 */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>

#ifdef __cplusplus
extern "C" {
#endif
#include "madx.h"
#ifdef __cplusplus
}
#endif

template<class T> // define class print operator
auto operator << (std::ostream& os, const T& t) -> decltype(t.Print(os), os)
{
  t.Print(os);
  return os;
}


namespace MaTh
{
  static int iMakeDipedge, iMakeEndMarkers, iMoreExpressions; // variables local to module that control makethin behavior
  static const int  el_type_maxlen=11;   // theElProp->get_el_type_maxlen();
  static const int  el_name_maxlen=25;
  static const int par_name_maxlen=19;
  static unsigned int Verbose;
  static const std::vector<std::string> DoNotCopy ={"l","lrad","slot_id","assembly_id","slice","comments"};
  static const std::vector<std::string> DoNotCopy2=           {"slot_id","assembly_id"};
  static const std::vector<std::string> WireCollimatorParmList={"xma","yma","closed_orbit", "current", "l_phy", "l_int"};
  static char ExtraChar='_';
}

//------------------------------- forward declarations --------------

class my_Element_List
{
public:
  my_Element_List(){}; // constructor
  ~my_Element_List()=default; // destructor
  std::vector<std::string> el_name_list;
  std::vector<element*> el_ptr;
  element* my_make_element(const std::string el_name, const std::string parent,command* def, int flag);
  void Print(std::ostream &StrOut = std::cout) const;
};

class ElmAttr // used for current element, list of all attributes defined and flag if on
{
public:
  ElmAttr(const element* el); // constructor
  ~ElmAttr()=default; // destructor
  void Print(std::ostream &StrOut = std::cout) const;
  void TurnOff(const std::vector<std::string>& off_list);
  void TurnOnActive(const element* el);
  std::vector<std::string> get_list_of_active_attributes() const;
  std::vector<std::string> AttrNam;
  std::vector<bool> On;
};

class SliceDistPos // defines the distances and positions, depending on number of slices and slicing style
{
public:
  SliceDistPos(const int n,const bool teapot_fl); // constructor
  ~SliceDistPos(); // destructor
  void Print(std::ostream &StrOut = std::cout) const;
  double delta;
  double Delta;
  //
  std::string delta_str,delta_half_str; // string expression
  std::string Delta_str,Delta_half_str; // string expression
private:
  int n; // number of slices
  bool teapot_fl;
};
SliceDistPos::~SliceDistPos()=default; // destructor

class OneElementWithSlices // One Element with Slices       used to work on slices, derived from thick_elem which is not modified - declared as constant
{
public:
  OneElementWithSlices(const element* thick_elem,element* sliced_elem); // constructor
  ~OneElementWithSlices()=default; // destructor
  const element* thick_elem; // pointer to the thick element
  std::vector<element*> theSlices; // pointer(s) to the one or several slices
};

class ElementListWithSlices
{
public:
  std::vector<OneElementWithSlices*> VecElemWithSlices; // vector of thick elements+slices
  ElementListWithSlices(unsigned int verbose); // constructor
  ~ElementListWithSlices(); // destructor
  void put_slice(const element* thick_elem, element* sliced_elem); // add sliced_elem to VecElemWithSlices
  element* find_slice(const element* thick_elem, const int slice); // find address of thin slice by slice number for thick_elem
  element* find_slice(const element* thick_elem, const std::string name); // find address of thin slice by slice name for thick_elem
  void Print(std::ostream &StrOut = std::cout) const;
  void PrintCounter(std::ostream &StrOut = std::cout) const;
private:
  int find_thick(const element* thick_elem); // find thick_element in VecElemWithSlices, <0 means not found,  used inside find_slice
  unsigned int verbose;
  unsigned int get_thin_calls,get_thin_iteractions; // to monitor the search (in-) efficiency find_slice
  int ilast1,ilast2; // keep last two found find_slice, useful in recursive searches which switch between slices and parents
};

class SequenceList
{
public:
  SequenceList(); // constructor
  ~SequenceList()=default; // destructor
  sequence* slice_sequence(const std::string slice_style,sequence* thick_sequ,const std::string LastSequenceSliced="",const std::string LastStyle="");
  void put_sequ(sequence* thick_sequ);       // add a sequence to the sequence list
  sequence* find_sequ(sequence* thick_sequ); // check if thick_sequ is already there, if yes return the pointer to it,  used to check if the sequence was already sliced
  void Print(std::ostream &StrOut = std::cout) const;
  void Reset();
private:
  std::vector<sequence*> my_sequ_list_vec; // list of sequences
};

static ElementListWithSlices *theSliceList=nullptr, *theRbendList=nullptr, *theBendEdgeList=nullptr; // global since MAD-X works with single global element_list
static my_Element_List *my_El_List=nullptr;

class SeqElList // sequence with elements considered for slicing
{
public:
  SeqElList(const std::string seqname,const std::string slice_style,/*sequence* thick_sequ,*/sequence* sliced_seq,node* thick_node,SequenceList* theSequenceList); // constructor
  ~SeqElList(); // destructor
  void Print(std::ostream &StrOut = std::cout) const;
  void slice_node(); // decides what to do : nothing, slice_node_translate, slice_node_default
  node* current_node() const { return thick_node;} // get
  void  current_node(node* thisnode) { work_node=thick_node=thisnode; } // set
private:
  double simple_at_shift(const int slices, const int slice_no) const;
  double teapot_at_shift(const int slices, const int slice_no) const;
  double collim_at_shift(const int slices, const int slice_no) const;
  double hybrid_at_shift(const int slices, const int slice_no) const;
  double at_shift(const int slices, const int slice_no,const std::string local_slice_style) const; // return at relative shifts from centre of unsliced magnet
  void kn_ks_from_thick_elem(const element* thick_elem,command_parameter* kn_pars[4],command_parameter* ks_pars[4]) const; // read k0-k3, k0s-k3s in thick_elem and put them in kn_pars, ks_pars
  void add_ktap(command_parameter* k_param,const element* thick_elem);
  void add_ktap_i(const int i,command_parameter* k_param,const std::string k_name,const std::string ktap_name,const element* thick_elem);
  command_parameter* make_k_list(const std::string parnam,command_parameter* k_pars[4]) const; // from k values 0-3 to  expr lists
  element*   new_marker_element(const std::string el_name, const element* el_inp);
  element*  create_wire_element(const element* thick_elem, int slice_no);
  element*        sbend_from_rbend(element* rbend_el);
  element*      create_thick_slice(const element* thick_elem,const int slice_type);
  element*      create_thin_slices(const element* thick_elem, int slice_no);
  element*    create_thin_solenoid(const element* thick_elem, int slice_no);
  element* create_thin_elseparator(const element* thick_elem, int slice_no);
  element*   create_sliced_element(const element* thick_elem, int slice_no);
  element* create_bend_dipedge_element(element* thick_elem,const bool Entry);
  void finish_make_sliced_elem(element*& sliced_elem, const element* thick_elem, command* cmd, const std::string parent_name, int slice_no); // final common steps
  void slice_node_translate();  // slice/translate and add slices to sliced sequence
  void slice_node_default();    // like collimator and add slices to sliced sequence
  void slice_attributes_to_slice(command* cmd,const element* thick_elem); // deal with attributes like kick
  void place_thick_slice(const element* thick_elem, element* sliced_elem, const int i);
  void place_start_or_end_marker(const bool at_start);
  node* copy_thin(node* thick_node);
  node* thick_node;  // current node, that is considered for slicing
  node*  work_node;  // clone of thick_node for non-destructive rbend->sbend translation
  const element* thick_elem_sliced; // last thick_elem sliced
  command_parameter *knl_param, *kns_param;
  SequenceList* theSequenceList;
  sequence *sliced_seq;
  std::string seqname; // name of the sequence
  std::string slice_style;
  unsigned int verbose;
  int nslices;
  const double eps;
  bool MakeDipedge; // translate dipoles   to    dipedge, dipole without edge effects, dipedge
};

//------------------------------- source code --------------

static const int        k_logical=0;
static const int            k_int=1;
static const int         k_double=2;
static const int        k_cstring=3;
static const int     k_int_array=11;
static const int  k_double_array=12;
static const int k_cstring_array=13;

static const bool dipedge_h1_h2_fl=false;   // normally false to avoid potentially non-simplectic partial higher order in dipedge. Optionally true as requested by Andrea Latina in 10/2014
static const bool kill_fringe_fl=true;      // requested by Laurent et al., somewhat redundant, should be sufficient to check existance of non-default h1,e1; h2,e2 parameters
static const bool Enable_all_attr_fl=true;  // set true to allow to enable all attibutes in the sliced sequence  --- otherwise only attributes defined in thick

// check general options
inline bool thin_foc_fl() { return get_option("thin_foc"); }
inline bool rbarc_fl()    { return get_option("rbarc"); } // by default on, then use (reduced) length of rbends

ElmAttr::ElmAttr(const element* el) // constructor
{
  if(el)
  {
    const command* el_cmd=el->def;
    const command_parameter_list *el_cmd_pl=el_cmd->par;
    for (int i=0;i<el_cmd_pl->curr; ++i)
    {
      AttrNam.push_back(el_cmd_pl->parameters[i]->name);
      On.push_back( el_cmd->par_names->inform[i] );
    }
    bool look_at_parent=true;
    element* el_parent= el->parent;
    while (look_at_parent && el_parent && el != el_parent && std::string(el_parent->name) != std::string(el_parent->base_type->name) )
    {
      TurnOnActive(el_parent);
      el_parent=el_parent->parent; // check recursively through parents
    }
  }
}

void ElmAttr::TurnOff(const std::vector<std::string>& off_list)
{
  if ( MaTh::Verbose>1 ) std::cout << "ElmAttr remove ";
  for(unsigned int i=0;i<AttrNam.size();++i)
  {
    std::string item=AttrNam[i];
    for(unsigned int j=0;j<off_list.size();++j)
    {
      if(item==off_list[j])
      {
        if ( MaTh::Verbose>1 ) std::cout << " " << item;
        On[i]=false;
        break;
      }
    }
  }
  if ( MaTh::Verbose>1 ) std::cout << std::endl;
}

void ElmAttr::TurnOnActive(const element* el)
{ // use when ElmAttr aleady existing for same size element, like parent
  if(el)
  {
    if ( MaTh::Verbose>1 ) std::cout << "ElmAttr turn on for " << el->name;
    const command* el_cmd=el->def;
    const command_parameter_list *el_cmd_pl=el_cmd->par;
    if( (int)AttrNam.size() == el_cmd_pl->curr)
    {
      for (int i=0;i<el_cmd_pl->curr; ++i)
      {
        if(el_cmd->par_names->inform[i] && !On[i])
        {
          On[i]=true;
          if ( MaTh::Verbose>1 ) std::cout << " " << AttrNam[i];
        }
      }
    }
    if ( MaTh::Verbose>1 ) std::cout << std::endl;
  }
}

void ElmAttr::Print(std::ostream &StrOut) const
{
  // show all attributes and mark in line below if on/off
  StrOut << std::right;
  StrOut << "ElmAttr "; for(unsigned int i=0;i<AttrNam.size();++i) StrOut << " " << AttrNam[i]; StrOut << '\n';
  StrOut << " On/off "; for(unsigned int i=0;i<AttrNam.size();++i) StrOut << " " << std::setw((int)AttrNam[i].length()) << On[i]; StrOut << std::endl;
}

std::vector<std::string> ElmAttr::get_list_of_active_attributes() const
{
  std::vector<std::string> result;
  for(unsigned int i=0;i<AttrNam.size();++i) if(On[i]) result.push_back(AttrNam[i]);
  return result;
}

static void warning_to_c(std::ostringstream& WarnStr) { warning((WarnStr.str()).c_str(),""); }

static bool NameIsInList(const std::string name,const std::vector<std::string>& namlist)
{
  for(int i=0; i < (int) namlist.size(); ++i)
    if( name==namlist[i] ) return true;
  return false;
}

static std::vector<std::string> str_v_join(const std::vector<std::string> v1,const std::vector<std::string> v2) // join two string vectors
{
  std::vector<std::string> result=v1;
  result.reserve(v1.size() + v2.size());
  for(unsigned int i=0;i<v2.size();++i) result.emplace_back(v2[i]);
  return result;
}

static double my_get_expression_value(expression* ex) // check for NULL and update the value as done in dump_expression
{
  if(ex)
  {
    //if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " ex->value=" << ex->value << " ex->status=" << ex->status << std::endl;
    ex->value = expression_value(ex, 2); // get/update value, makes sure the value stored agrees with the expression
    //if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " ex->value=" << ex->value << " ex->status=" << ex->status << std::endl;
  }
  return ex->value; // no expression, just value
}

static std::string my_dump_expression(expression* ex,bool detailed=false) // dump_expression in mad_expr.c only prints the value, here show also the expression and check for NULL
{
  std::ostringstream ostr;
  // if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << std::endl;
  ostr << std::setprecision(15) << "expression ";
  if(ex==nullptr) ostr << " is nullptr";
  else
  {
    if(ex->string) ostr << " string=" << std::left << std::setw(MaTh::el_name_maxlen) << ex->string << std::right;
    if(detailed) ostr << " status=" << ex->status;
    if(detailed) ostr << " polish=" << ex->polish;
    if(detailed) ostr << " value=" << std::setw(8) << ex->value;
    if(detailed) ostr << " my_get_expression_value=" << std::setw(8) << my_get_expression_value(ex);
    else ostr << " value=" << std::setw(8) << my_get_expression_value(ex);
  }
  return ostr.str();
}

static expression* expr_from_value(double value)
{
  std::ostringstream ostr;
  ostr << std::setprecision(15) << value; // use the value as string
  char exprstring[80];
  strcpy(exprstring,ostr.str().c_str());
  expression* expr = new_expression(exprstring, deco); // where deco is a global. Seems expression value not well defined
  // expression* expr = new_expression(exprstring, nullptr); // where deco is a global. Seems expression value not well defined
  // expression* expr = new_expression(ostr.str().c_str(), nullptr); // where deco is a global. Seems expression value not well defined
  expr->value=value;
  expr->status=1;
  return expr;
}

static expression* expr_from_value_2(double value)
{
  std::ostringstream ostr;
  ostr << std::setprecision(15) << value; // use the value as string
  expression* new_expr = new_expression(ostr.str().c_str(), deco); // where deco is a global.
  new_expr = compound_expr( new_expr,0,"+",0,0,1); // trick to update value
  return new_expr;
}

static double my_get_int_or_double_value(const element* el,const std::string parnam,bool &found) // works for integer and double, also useful as   my_get_int_or_double_value(el->base_type,char* parnam);  to get the default values
{
  // just returning  el_par_value(parnam,base_el);   or  el_par_value_recurse(parnam,base_el);
  //      is not good enough, gets 0 for integer parameters and in some cases too specific - like checks for non-zero length of dipoles  which does not allow to compare with the base type for dipoles
  // rather do all directly here, descending from el / el->def / el->def->par to parameters[i], loop through them and look at integer and double values
  // in case of expression uses the expression_value
  double val=0;
  found=false;
  if(el && el->def && el->def->par)
  {
    command_parameter_list* pl=el->def->par;
    for (int i = 0; i < pl->curr; ++i)
    {
      if(pl->parameters[i])
      {
        command_parameter* cp=pl->parameters[i];
        if( std::string(cp->name) == parnam )
        {
          if(cp->expr)
          {
            val = my_get_expression_value(cp->expr);
            found=true;
          }
          else switch (cp->type)
          {
            case k_int:    //    int value of expression, actually same as double value
              found=true;
              val = cp->double_value;
              break;
            case k_double: // double value of expression
              found=true;
              val = cp->double_value;
              break;
          }
        }
      }
    }
  }
  return val;
}

static double cmd_par_val(const command_parameter* par) // return the double value given by expression or directly the value
{
  double result=0;
  if(par)
  {
    if(par->type==k_double)
    {
      if(par->expr) result=expression_value(par->expr,0); else result=par->double_value;
    }
  }
  return result;
}

static std::string my_dump_command_parameter(const command_parameter* cp) // dump_command_parameter in mad_cmdpar.c only prints the value, here show also the expression by calling my_dump_expression and check for NULL
{
  std::ostringstream ostr;
  ostr << std::setprecision(15) << "my_dump_command_parameter ";
  if(cp==nullptr) ostr << " is nullptr" << '\n';
  else
  {
    ostr << "parameter:" << std::left << std::setw(MaTh::par_name_maxlen) << cp->name;
    ostr << std::right << " cp->type=" << std::setw(2) << cp->type;
    ostr << " stamp=" << cp->stamp << " ";
    double default_val=0;
    const double eps=1.e-15; // used to check if strength is compatible with zero
    switch (cp->type)
    {
      case k_logical:
        ostr << "logical: ";
        if( (int) cp->double_value) ostr << "true"; else ostr << "false";
        ostr << std::endl;
        break;
      case k_int:    // int    value of expression, actually same as double value
      case k_double: // double value of expression
        if(cp->expr) ostr << my_dump_expression(cp->expr); else ostr << " expression=nullptr ";
        if(cp->call_def) default_val=cp->call_def->double_value;
        if(cp->expr==nullptr && fabs(cp->double_value-default_val)>eps)
          ostr << " value=" << std::setw(10) << cp->double_value << " default=" << std::setw(10) << default_val;
        ostr << std::endl;
        break;
      case k_int_array:    // int array,     expr_list
      case k_double_array: // double array,  expr_list, used for example for Aperture, http://mad.web.cern.ch/mad/madx.old/Introduction/aperture.html
        if (cp->double_array)
        {
          if (cp->expr_list) // calculate the values
          {
            ostr << "array of " << cp->double_array->curr << "  ";
            for (int ei = 0; ei < cp->double_array->curr; ++ei)
            {
              if (ei < cp->expr_list->curr && cp->expr_list->list[ei] != nullptr)
              {
                ostr << std::right << std::setw(3) << ei << " :" << std::left << my_dump_expression(cp->expr_list->list[ei]) << std::right; // show expression and value
              }
              else ostr << std::right << std::setw(3) << ei << " is nullptr ";
            }
            ostr << "double array:";
            for (int ei = 0; ei < cp->double_array->curr; ++ei) ostr << cp->double_array->a[ei] << " ";
          }
        }
        ostr << std::endl;
        break;
      case k_cstring:
        ostr << "cstring:";
        if(cp->string) ostr << cp->string; else ostr << " NULL";
        ostr << std::endl;
        break;
      case k_cstring_array: // string array
        dump_char_p_array(cp->m_string);
        /* FALLTHRU */

      case '?':
        ostr << " cp->type=" << cp->type << " no info dump implemented so far" << '\n';
    }
  }
  return ostr.str();
}

static std::string my_dump_command_parameter_list(command_parameter_list* pl)
{
  std::ostringstream ostr;
  ostr << std::setprecision(15) << "my_dump_command_parameter_list";
  if(pl)
  {
    ostr << " name=" << std::setw(MaTh::par_name_maxlen) << pl->name;
    ostr << " curr=" << pl->curr << " max=" << pl->max << '\n';
    if(pl->curr > pl->max)
    {
      ostr << "*** error *** seen in my_dump_command_parameter_list max=" << pl->curr << " > " << " curr" << pl->curr << " set curr back to max" << '\n';
      pl->curr = pl->max;
    }
    for (int i = 0; i < pl->curr; ++i)
    {
      ostr << std::setw(2) << i << " : ";
      if(pl->parameters[i]) ostr << my_dump_command_parameter(pl->parameters[i]); else ostr << " NULL ";
    }
  }
  else ostr << " is NULL";
  return ostr.str();
}

static void SetParameterValue(const std::string parnam,element* el,const double val,const int type=k_double) // set value and type, by default double
{
  command*   el_def=el->def;
  name_list* nl=el_def->par_names;
  const int ei=name_list_pos(parnam.c_str(),nl);
  if(ei > -1)
  {
    command_parameter* cp=el_def->par->parameters[ei];
    if(cp)
    {
      if(MaTh::Verbose>1)
      {
        std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " el->name=" << std::setw(MaTh::el_name_maxlen) << el->name << " parameter " << parnam
        << " was double_value=" << cp->double_value
        << " and type=" << cp->type;
        if(cp->expr) std::cout << " has " << my_dump_expression(cp->expr); else std::cout << " no expression";
        std::cout << " set to val=" << val
        << " and type=" << type << '\n';
      }
      if(cp->expr) cp->expr=nullptr; // remove any expression
      cp->double_value=val; // set the double value
      cp->type=type; // set the type value
    }
  }
  else
  {
    std::ostringstream WarnStr;
    WarnStr << "SetParameterValue for parameter " << parnam << " failed for " << std::setw(MaTh::el_name_maxlen) << el->name << " parameter not in element name_list";
    warning_to_c(WarnStr);
  }
}

static void SetParameter_in_cmd(command* cmd, const command_parameter* param, const std::string parnam,const int inf)
{ // set param in cmd if exists
  // if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " parnam=" << parnam << " param=" <<  param << '\n';
  if(param)
  {
    name_list* nl=cmd->par_names;
    const int ei=name_list_pos(parnam.c_str(),nl);
    if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << std::left << " parnam=" << std::setw(MaTh::par_name_maxlen) << parnam << " param=" <<  param << " parame->name=" << std::setw(MaTh::par_name_maxlen) << param->name << " ei=" << std::setw(2) <<  ei << std::right << '\n';
    if(ei > -1)
    {
      nl->inform[ei]=inf;
      command_parameter* param_copy=clone_command_parameter(param); // make copy which can be modified
      if( parnam!=std::string(param->name) ) strcpy(param_copy->name, parnam.c_str()); // use parnam, can be different from parnam->name,  like e2 of bend which becomes e1 in exit dipedge
      cmd->par->parameters[ei]=param_copy;
    }
  }
}

static void ParameterTurnOn(const std::string parnam,element* el) // request that this parameter is written to output
{
  const command* el_def=el->def;
  name_list* nl=el_def->par_names;
  const int ei=name_list_pos(parnam.c_str(),nl);
  if(ei > -1) nl->inform[ei]=1; // Turn on by setting inform to 1
  else
  {
    std::ostringstream WarnStr;
    WarnStr << "ParameterTurnOn for parameter " << parnam << " failed for " << std::setw(MaTh::el_name_maxlen) << el->name << " parameter not in element name_list ";
    warning_to_c(WarnStr);
  }
}

static void ParameterRemove(const std::string parnam,element* el)
{
  if(el)
  {
    const command* el_def=el->def;
    name_list* nl=el_def->par_names;
    const int ei=name_list_pos(parnam.c_str(),nl);
    if(ei > -1)
    {
      nl->inform[ei]=0; // Turn off  by setting inform to 0   -- effective in save,  but still used in twiss; so really delete expression and turn value off
      double default_value=0;
      default_value=el->base_type->def->par->parameters[ei]->double_value;  // el_par_value(parnam,el->base_type) cannot be used here, base element length may be zero
      command_parameter* cp=el_def->par->parameters[ei];
      if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " in " << el-> name << " parameter" << std::setw(12) << parnam
        << " value=" << std::setw(6) << cp->double_value << " set to default=" << std::setw(6) << default_value
        << " for " << std::setw(12) << parnam << " cp->expr=" << cp->expr << " and set expression to nullptr" << '\n';
      cp->type = k_double;
      cp->double_value = default_value;
      cp->expr=nullptr;
    }
  }
}

static command* new_cmdptr(const element* elem)
{ // new command* with all attributes as defined in elem
  command* cmd=clone_command(elem->def); // definition has all inform on
  for(int i=0; i< cmd->par->curr; ++i) cmd->par_names->inform[i]=0; // turn off
  return cmd;
}

static std::string my_dump_command(const command* cmd)
{
  std::ostringstream ostr;
  if(cmd==nullptr) ostr << " is NULL";
  else
  { // command defined in mad_elem.h, c-structures based on pointing to pointers, contains name_list par_names and command_parameter_list par
    ostr << "my_dump_command command: ";
    ostr << cmd->name;
    ostr << "  module: "   ; ostr << cmd->module;
    ostr << "  group: "    ; ostr << cmd->group;
    ostr << "  stamp= " << cmd->stamp << "  link_type= " << cmd->link_type << "  mad8_type= " << cmd->mad8_type;
    ostr << "  #par_names "; if(cmd->par_names->curr)   ostr << cmd->par_names->curr; else ostr << " NULL";
    ostr << "  #par= "     ; if(cmd->par->curr)   ostr << cmd->group; else ostr << " NULL";
    ostr << '\n';
    for (int i = 0; i < cmd->par->curr; ++i) ostr << my_dump_command_parameter(cmd->par->parameters[i]);
    ostr << "within command par:";       if(cmd->par)       ostr << '\n' << my_dump_command_parameter_list(cmd->par); else ostr << " NULL" << '\n';
  }
  ostr << '\n';
  ostr << "my_dump_command command end" << std::endl;
  return ostr.str();
}

static std::string my_dump_element(const element* el)
{
  std::ostringstream ostr;
  ostr << std::setprecision(15) << std::left << '\n' << "my_dump_element";
  if(el==nullptr) ostr << " is NULL";
  else
  { // element defined by c-structures based on pointing to pointers
    ostr << " name=" << std::setw(MaTh::el_name_maxlen) << el->name;
    if(el->base_type) ostr << " base_type=" << el->base_type->name;
    ostr << " stamp=" << el->stamp << " length=" << el->length << " parent name=" << std::setw(MaTh::el_type_maxlen) << el->parent->name;
    ostr << " def_type=" << el->def_type;
    if(el->def_type) ostr << " which means defined separately"; else ostr << " which means inside sequence";
    ostr << '\n';
    ostr << "within element " << my_dump_command(el->def);
  }
  return ostr.str();
}

static void Remove_All_Fringe_Field_Parameters(element* el)
{
  static std::vector<std::string> FringePar = {"e1","e2","fint","fintx","h1","h2","hgap"};
  for(unsigned int i=0; i < FringePar.size(); ++i) ParameterRemove(FringePar[i].c_str(),el);
  if(kill_fringe_fl)
  {
    SetParameterValue("kill_ent_fringe",el,true,k_logical);
    SetParameterValue("kill_exi_fringe",el,true,k_logical);
    ParameterTurnOn("kill_ent_fringe",el); // turn writing on
    ParameterTurnOn("kill_exi_fringe",el); // turn writing on
  }
  if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " el name=" << std::left << std::setw(MaTh::el_name_maxlen) << el->name << " base_type" << el->base_type->name << " after  remove : " << my_dump_element(el) << '\n';
}

static std::string my_dump_node(const node* node)
{
  std::ostringstream ostr;
  ostr << std::setprecision(15) << '\n' << "my_dump_node";
  if(node==nullptr) ostr << " is NULL";
  else
  {
    ostr << std::setprecision(15) << " node:";
    std::string pname, nname, from_name;
    if (node->previous  != nullptr)     pname = node->previous->name;
    if (node->next      != nullptr)     nname = node->next->name;
    if (node->from_name != nullptr) from_name = node->from_name;
    ostr << " name=" << std::left << std::setw(MaTh::el_name_maxlen) << node->name << std::right
    << " occ=" << node->occ_cnt
    << " node->base_name=" << std::left << std::setw(MaTh::el_type_maxlen) << node->base_name << std::right
    << " from_name=" << std::left << std::setw(10) << from_name << std::right
    << " at_value=" << std::setw(10) << node->at_value
    << " position=" << std::setw(10) << node->position
    << " previous=" << std::setw(MaTh::el_name_maxlen) << pname
    << " next=" << std::setw(MaTh::el_name_maxlen) << nname << std::right
    << " at_expr: ";
    if(node->at_expr) ostr << my_dump_expression(node->at_expr); else ostr << "NULL ";
    if(node->p_elem) ostr << my_dump_element(node->p_elem);
    if(node->cl) { for(int i=0; i< node->cl->curr; ++i) dump_constraint(node->cl->constraints[i]); }
  }
  ostr << std::endl;
  return ostr.str();
}

static std::string my_dump_sequence(const sequence* c_sequ,const int level)
{ // level 1 little info, 3 dump also nodes, 4 dump also elements
  std::ostringstream ostr;
  if(c_sequ==nullptr) ostr << "sequence is NULL";
  else
  {
    node* c_node;
    double suml = zero;
    ostr << "sequence:" << c_sequ->name;
    if(c_sequ->refpos)    ostr << " refpos=" << c_sequ->refpos; else ostr << " refpos=NULL";
    if(c_sequ->next_sequ) ostr << " next_sequ=" << c_sequ->next_sequ; else ostr << " next_sequ=NULL";
    ostr << " ref_flag=" << c_sequ->ref_flag;  // -1 for exit, 0 for centre, 1 for entry
    if(c_sequ->ref_flag==-1) ostr << " (exit) ";
    else if(c_sequ->ref_flag==0) ostr << " (centre) ";
    else if(c_sequ->ref_flag==1) ostr << " (entry) ";
    ostr << " share=" << c_sequ->share << " nested=" << c_sequ->nested << " con_cnt=" << c_sequ->con_cnt << " stamp=" << c_sequ->stamp << " line=" << c_sequ->line << " add_pass=" << c_sequ->add_pass << " length=" << c_sequ->length << '\n';
    c_node = c_sequ->start;
    ostr << std::setprecision(15);
    double lastvalue=0;
    while(c_node != nullptr)
    {
      suml += c_node->length;
      if (level > 2)
      {
        ostr << my_dump_node(c_node);
        if (level > 3 && c_node->p_elem != nullptr) ostr << my_dump_element(c_node->p_elem);
      }
      else if (level > 0 && strcmp(c_node->base_name, "drift") != 0)
      {
        ostr << std::left << std::setw(MaTh::par_name_maxlen) << c_node->name << std::right
        << " at_value=" <<  std::setw(10) << c_node->at_value
        << " position=" <<  std::setw(6) << c_node->position
        << " length="   << std::setw(17) << c_node->length;
        if(c_node->p_elem) // show lrad if available
        {
          double lrad=el_par_value("lrad",c_node->p_elem);
          if(lrad>0) ostr << " lrad value=" << lrad;
        }
        if(c_node->from_name) ostr << " from " << c_node->from_name;
        if(c_node->at_expr) ostr << " at_expr " << my_dump_expression(c_node->at_expr);
        if(c_node->at_expr) { double currentvalue=my_get_expression_value(c_node->at_expr); ostr << " diff=" << std::setw(8) << currentvalue-lastvalue; lastvalue=currentvalue; }
        ostr << '\n';
      }
      if (c_node == c_sequ->end)  break;
      c_node = c_node->next;
    }
    ostr << "===== sum of node length=" << std::setw(8) << suml << '\n';
    ostr << std::endl;
  }
  return ostr.str();
}

static std::string my_get_cmd_expr_str(const command_parameter* cmd) // return the expression as string, if there is only a value, return the value as string
{
  std::string result="";
  if(cmd)
  {
    if(cmd->expr && cmd->expr->string) result=cmd->expr->string; // the expression is define as string, use this
    else // look for a value
    {
      const double eps=1.e-15; // used to check if strength is compatible with zero
      if( fabs(cmd->double_value)>eps ) // value defined and non-zero
      {
        std::ostringstream ostr;
        if(cmd->double_value<0) ostr << "("; // enclose negative value in brackets
        ostr << cmd->double_value;
        if(cmd->double_value<0) ostr << ")"; // enclose negative value in brackets
        result=ostr.str();
      }
    }
  }
  if(result.length()<1) result="0";
  if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " result=" << result << '\n';
  return result;
}

static expression* my_get_param_expression(const element* el,const std::string parnam) // get a new copy of the expression for a parameter from an element, use the value as new expression if the expression was NULL
{
  const command* el_def=el->def;
  name_list* nl=el_def->par_names;
  const int ei=name_list_pos(parnam.c_str(),nl);
  expression* ep=nullptr;
  if(ei > -1)
  {
    const command_parameter* cp = el_def->par->parameters[ei]; // pointer to original parameter
    if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " for element " << std::setw(MaTh::par_name_maxlen) << el->name << " parameter " << std::setw(MaTh::par_name_maxlen) << parnam << " ei=" << ei << " my_dump_expression(cp->expr):" << my_dump_expression(cp->expr) << " cp->double_value=" << cp->double_value << '\n';
    command_parameter* cp_copy = clone_command_parameter( cp );  // copy of the original parameter that can be modified
    if(cp_copy->expr==nullptr) cp_copy->expr = expr_from_value(cp->double_value);
    else if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " expression exists use it, expr=" << my_dump_expression(cp_copy->expr) << '\n';
    ep = cp_copy->expr;
  }
  return ep;
}

static expression* my_get_param_expression(const command_parameter* cp) // get a new copy of the expression for a parameter, use the value as new expression if the expression was NULL
{
  expression* ep=cp->expr;
  if(ep==nullptr) // no expression yet
  { // use the value as new expression if the expression was NULL
    std::ostringstream ostr;
    ostr << std::setprecision(15) << cp->double_value; // use the value as string
    ep = new_expression(ostr.str().c_str(),deco); // where deco is a global.  // cp_copy->expr->value = cp->double_value; // seems to have no effect and this not needed
    if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " create new expression from string " << ostr.str() << " now " << my_dump_expression(ep) << '\n';
  }
  return ep;
}

element* my_Element_List::my_make_element(const std::string el_name, const std::string parent,command* def, int flag)
{
  for(unsigned int i=0;i<el_name_list.size();++i)
  { // check if already there, not done in MAD-X make_element
    if( el_name_list[i]==el_name)
    {
      element* elm=el_ptr[i]; // candidate, look at base
      if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " el_name=" << std::setw(10) << el_name << " base=" << elm->base_type->name << " base wanted=" << std::setw(10) << parent<< " already done " << '\n';
      return elm;
    }
  }
  element* new_elm=make_element(el_name.c_str(),parent.c_str(),def,flag);
  el_name_list.push_back(el_name);
  el_ptr.push_back(new_elm);
  return new_elm;
}

void my_Element_List::Print(std::ostream &StrOut) const
{
  StrOut << std::right;
  for(unsigned int i=0;i<el_name_list.size();++i)
  {
    StrOut << std::setw(5) << i << std::setw(40) << el_name_list[i] << " " << el_ptr[i] << std::endl;
  }
}

static command_parameter* par_scaled(const command_parameter* par_inp, const command_parameter* length_param, const std::string new_par_name, const int nslices)
{ // scale parmeter  * length / nslices      and give new name
  command_parameter* par_out=nullptr;

  if( par_inp && length_param )
  {
    par_out=clone_command_parameter(par_inp); // start from clone of input parameter
    strcpy(par_out->name,new_par_name.c_str());
    if (length_param->expr ) // first * l, expression or value
    {
      if(!par_out->expr) par_out->expr=my_get_param_expression(par_out);
      par_out->expr = compound_expr(par_out->expr,par_out->double_value,"*",length_param->expr,length_param->double_value,1); // multiply expression with length
    }
    else
    {
      // if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " multiply par_out->double_value=" << par_out->double_value << '\n';
      par_out->double_value *= length_param->double_value; // multiply value with length
      // if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << "      now par_out->double_value=" << par_out->double_value << '\n';
    }
    if (nslices > 1) // 2nd step, divide by nslices, expression or number
    {
      if (par_out->expr) par_out->expr = compound_expr(par_out->expr,0.,"/",nullptr,nslices,1);
      else par_out->double_value /= nslices;
    }
    if( MaTh::iMoreExpressions<1 && par_out->expr)
    {
      par_out->double_value=my_get_expression_value(par_out->expr);
      par_out->expr=nullptr;
    }
   }
  return par_out;
}

static bool thick_fl(const element* el) // true if the element has a thick parameter and if the value is positive, 0 otherwise
{
  const command* el_def=el->def;
  name_list* nl=el_def->par_names;
  const int ei=name_list_pos("thick",nl);
  return (ei > -1 && el_def->par->parameters[ei]->double_value > 0);
}

static std::string dump_slices(el_list* the_element_list) // Loops over all current elements and prints the number of slices. Used for debug and info
{
  std::ostringstream ostr;
  ostr << "++++++ dump_slices";
  if(MaTh::Verbose>1)
    ostr << "   verbose on, all elements are listed" << '\n';
  else
    ostr << "   only elements with non default selection (other than 1 thin) are shown" << '\n';
  ostr << "            name  #slices      derived from  #slices" << '\n';
  int n_elem_with_slice=0,n_elem_with_slice_gt_1=0;
  for(int i=0; i< the_element_list->curr; ++i) // loop over the_element_list
  {
    element* el = the_element_list->elem[i];
    const command* el_def=el->def;
    name_list* nl=el_def->par_names;
    const int ei=name_list_pos("slice",nl);
    if(ei > -1) // element with slice number defined
    {
      n_elem_with_slice++;
      const int slices=el->def->par->parameters[ei]->double_value;
      int slices_parent=0;
      std::string parent_name="no parent";
      if(el->parent!=nullptr) // look also at parent if existing
      {
        slices_parent=el->parent->def->par->parameters[ei]->double_value;
        parent_name=el->parent->name;
      }
      if(slices>1) n_elem_with_slice_gt_1++;
      if(MaTh::Verbose>1|| slices !=1 || thick_fl(el)) // print all with verbose. with debug skip elements with default selection thin 1 slice
      {
        ostr << " " << std::setw(15) << el->name << " " << slices;
        if(thick_fl(el)) ostr << " thick"; else ostr << " thin ";
        if(el != el->parent)
        {
          ostr << " " << std::setw(18) << parent_name << " " << slices_parent; // show also parent if not same as child
          if(thick_fl(el->parent)) ostr << " thick"; else ostr << " thin ";
        }
        ostr << std::endl;
      }
    }
  }
  if(MaTh::Verbose>1) ostr << "       general option thin_foc=" << thin_foc_fl() << '\n'; // global option like debug or verbose, not element specific, still print here for info
  ostr << "------ end of dump slices. There were " << std::setw(4) << the_element_list->curr << " elements, "
  << std::setw(3) << n_elem_with_slice << " with slice numbers and "
  << std::setw(2) << n_elem_with_slice_gt_1 << " with slice numbers>1\n\n";
  return ostr.str();
}

static void force_consistent_slices(el_list* the_element_list) // hbu 10/2005 loop over all elements and check that #slices of child and parent agree,  if not, use the maximum for both
{
  for(int i=0; i< the_element_list->curr; ++i) // loop over the_element_list
  {
    element* el = the_element_list->elem[i];
    const command* el_def=el->def;
    name_list* nl=el_def->par_names;
    const int ei=name_list_pos("slice",nl);
    if(ei > -1 && el->parent!=nullptr && el != el->parent )
    {
      command_parameter*  child=el_def->par->parameters[ei];
      command_parameter* parent=el->parent->def->par->parameters[ei];
      int slices=child->double_value;
      int slices_parent=parent->double_value;
      if(slices != slices_parent)
      {
        if(slices>slices_parent) slices_parent=slices; else slices=slices_parent;
        child->double_value=slices;
        parent->double_value=slices_parent;
        int el_thick_pos = name_list_pos("thick",nl);
        if(el_thick_pos > -1) el->parent->def->par->parameters[el_thick_pos]->double_value = el_def->par->parameters[el_thick_pos]->double_value; // copy thick flag from child to parent
      }
    }
  }
}

static int get_slices_from_elem(const element* el)
{
  int slices=1; // default
  if(el)
  {
    name_list* nl=el->def->par_names;
    const int ei=name_list_pos("slice",nl);
    if( ei > -1) slices=el->def->par->parameters[ei]->double_value;
    if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__  << " element " << el->name << " ei=" << ei << " slices=" << slices << '\n';
  }
  return slices;
}

static std::string make_thin_name(const std::string name_inp,const int slice) // make a node name from element name and slice number
{ // example     e_name=mqxa.1r1 slice=1 result=mqxa.1r1..1
  std::string name_out=name_inp+".."+std::to_string(slice);
  if(name_out.length()>NAME_L)
  {
    warning(std::string("slice name is too long, truncated at "+std::to_string(NAME_L)+" characters").c_str(), name_out.c_str());
    name_out=name_out.substr(0,NAME_L);
  }
  return name_out;
}

static command_parameter*
scale_and_slice(command_parameter* knsl_param,const command_parameter* length_param,const int nslices,const bool mult_with_length) // multiply the kn or ks parameter by the length and divide by slice
{
  int last_non_zero=-1;
  if (knsl_param   == nullptr) return nullptr;
  if (length_param == nullptr) return nullptr;
  const double eps=1.e-15; // used to check if strength is compatible with zero

  for (int i=0; i<knsl_param->expr_list->curr; ++i)
  {
    expression* kn_i_expr = knsl_param->expr_list->list[i];
    double kn_i_val  = knsl_param->double_array->a[i];

    if ((kn_i_expr != nullptr && zero_string(kn_i_expr->string)==0) || fabs(kn_i_val)>eps )
    {
      last_non_zero=i;
      if ( mult_with_length||i>0 ) // apply mult_with_length only to zero order multipole
      {
        if (length_param->expr || kn_i_expr)
          kn_i_expr = compound_expr(kn_i_expr, kn_i_val, "*", length_param->expr, length_param->double_value,1); // multiply expression with length
        else kn_i_val *= length_param->double_value; // multiply value with length
      }
      if (nslices > 1) // give the correct weight by slice (divide by the number of slices)
      {
        if (kn_i_expr) kn_i_expr = compound_expr(kn_i_expr,kn_i_val,"/",nullptr,nslices,1);
        else kn_i_val *= 1./nslices;
      }
    }
    if(kn_i_expr) knsl_param->expr_list->list[i] = kn_i_expr;
    knsl_param->double_array->a[i] = kn_i_val;
  } // for i ..
  if (last_non_zero==-1)
  {
    delete_command_parameter(knsl_param);
    knsl_param=nullptr;
  }
  if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " end" << std::endl;
  return knsl_param;
}

static void set_lrad(command* cmd,const command_parameter* length_param,const int nslices)
{ // to keep information of original length divided by nubmer of nslices
  if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " length_param=" << length_param << " nslices=" << nslices  << '\n';
  if(length_param)
  {
    name_list* nl=cmd->par_names;
    const int ei=name_list_pos("lrad",nl);
    if(ei>-1)
    {
      if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " for lrad  ei=" <<  ei << '\n';
      command_parameter* l_par = cmd->par->parameters[ei] = clone_command_parameter(length_param); // keep what was l
      strcpy(l_par->name,"lrad"); // rename to lrad
      if (nslices > 1) // and divide numbers or expressions by the number of nslices
      {
        if (l_par->expr) l_par->expr = compound_expr(l_par->expr,0.,"/",nullptr,nslices,1);
        else l_par->double_value /= nslices;
      }
      if( MaTh::iMoreExpressions<1 && l_par->expr)
      {
        double val=my_get_expression_value(l_par->expr);
        l_par->double_value=val;
        l_par->expr=nullptr;
      }
      nl->inform[ei]=1;
    }
    else if (MaTh::Verbose) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " *** warning *** , element has no lrad, nothing done" << " length_param=" << length_param << " nslices=" << nslices  << '\n';
  }
}

static expression* curved_from_straight_length(const element* rbend_el)
{
  expression* l_rbend_expr = my_get_param_expression(rbend_el,"l"); // get expression or create new from constant
  expression* l_sbend_expr = nullptr;
  if(rbarc_fl() && l_rbend_expr ) // increase the straight rbend length to sbend length
  { // RBEND, "l" parameter  el_par_value("l","rbend") with rbarc on as default is the shorter straight length,   val = l * angle / (two * sin(angle/two));
    // this is also shown in twiss node and element length give always the curved length
    // in going from rbend to sbend, this correction must be applied if the "l" expression is used, not for the value
    std::string anglestr = my_get_cmd_expr_str( return_param_recurse("angle", rbend_el) );
    const std::string rat = "1.0/sinc("+anglestr+"*0.5)"; // L_sbend / L_rbend
    expression* rat_expr = new_expression(rat.c_str(),deco);
    l_sbend_expr = compound_expr(l_rbend_expr,0,"*",rat_expr,0,1); // this also updates the value
    if(MaTh::Verbose>1)
    {
      bool found=false;
      double straight_length=my_get_int_or_double_value(rbend_el,"l",found);
      std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__  << " " << rbend_el->name << " rbarc on, increase rbend straight length expression of value " << straight_length << " to curved sbend length  using anglestr=" << anglestr
      << " updated l_sbend_expr " << my_dump_expression(l_sbend_expr) << " value should now be same as the curved rbend_el->length=" << rbend_el->length << '\n';
    }
  }
  else l_sbend_expr=l_rbend_expr;
  return l_sbend_expr;
}

static void add_node_at_end_of_sequence(node* node,sequence* sequ) // position in thin sequence defined with at_value and from_name
{
  if (sequ->start == nullptr) // first node in new sequence
  {
    sequ->start = node;
    node->next = nullptr;
    node->previous = nullptr;
  }
  else // add to end
  {
    sequ->end->next = node;
    node->previous  = sequ->end;
  }
  sequ->end = node;

  if(MaTh::Verbose>1)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << std::left << std::setw(MaTh::el_name_maxlen) << node->name << " " << std::setw(MaTh::par_name_maxlen) << node->base_name << std::right
    << " position=" << std::setw(10) << node->position  << " at_value=" << std::setw(10) << node->at_value;
    if(node->at_expr)   std::cout << " " << my_dump_expression(node->at_expr);
    if(node->from_name) std::cout << " from " << std::setw(5) << node->from_name; else std::cout << "           ";
    std::cout << " length="   << std::setw(10) << node->length  << " in " << sequ->name << '\n';
  }
  add_to_node_list(node, 0, sequ->nodes);
  return;
}

static void add_half_angle_to(const element* rbend_el,element* to_el,const std::string to_parm_name) // get half surface angle of rbend, add to to_parm_name (e1 or e2) of dipedge or sbend
{
  if(rbend_el && to_el)
  {
    expression* half_angle_expr  = compound_expr(my_get_param_expression(rbend_el,"angle"),0.,"/",nullptr,2,1); // angle/0.5  add this to any existing to_parm_name (e1 or e2)
    command_parameter* to_param = return_param_recurse(to_parm_name.c_str(),to_el); // get param from element, use here the non const version of return_param_recurse, to modify to_el
    if(to_param) // modify the existing parameter in to_el
    {
      if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " to_parm_name=" << to_parm_name << "    original to_param " << my_dump_expression(my_get_param_expression(to_param)) << '\n';
      to_param->expr = compound_expr( my_get_param_expression(to_param) ,0,"+",half_angle_expr,0,1);
    }
    else // param was null in to_el, start from parameter definition
    {
      int ipar = name_list_pos(to_parm_name.c_str(),to_el->def->par_names);
      if(ipar > -1) // already in name_list
      {
        to_param = clone_command_parameter( to_el->def->par->parameters[ipar] );  // copy of the original length parameter that can be modified
        to_el->def->par->parameters[ipar]->expr=half_angle_expr;
        ParameterTurnOn(to_parm_name,to_el);
        if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " use existing to_param from ipar= " << ipar
          << " to_param=" << to_param
          << " to_el->def->par->parameters[ipar]->expr=" << to_el->def->par->parameters[ipar]->expr << '\n';
      }
      to_param->expr = half_angle_expr;
    }
    if( MaTh::iMoreExpressions<1 && to_param->expr)
    {
      to_param->double_value=my_get_expression_value(to_param->expr);
      to_param->expr=nullptr;
    }
  }
}

static void place_node_at(const node* node, sequence* to_sequ, element* sliced_elem,expression* at_expr)
{
  struct node* this_node = new_elem_node(sliced_elem, node->occ_cnt);
  double at = node->at_value;
  this_node->from_name = node->from_name;
  this_node->at_value  = at;
  if(at_expr) this_node->at_expr   = at_expr;
  if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " place " << sliced_elem->name << " using at_expr where " << my_dump_expression(at_expr) << " at_value=" << at << std::endl;
  add_node_at_end_of_sequence(this_node,to_sequ); // add the thick node to the sequences
  if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " dump_node(this_node)   :" << std::endl;
}

static void place_thin_slice(const node* node, sequence* to_sequ, element* sliced_elem,const double rel_shift) // used to place dipedge and markers
{
  if(node->p_elem)
  {
    if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " MaTh::iMoreExpressions=" << MaTh::iMoreExpressions
      << " sliced_elem " << sliced_elem->name << '\n';
    double at = node->at_value;
    expression* length_param_expr=my_get_param_expression(node->p_elem, "l"); // get expression or create new from constant
    expression* at_expr=nullptr;
    if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " sliced_elem=" << sliced_elem->name << " node->p_elem=" << node->p_elem->name << " length_param_expr " << my_dump_expression(length_param_expr) << " node->at_expr " << my_dump_expression(node->at_expr) << " rel_shift=" << rel_shift << '\n';
    if( MaTh::iMoreExpressions<1 )
    {
      // at_expr = compound_expr(node->at_expr, at, "+",  nullptr, my_get_expression_value(length_param_expr) *rel_shift );  // use length and shift values, no expressions; failed to update the value
      double at_shift=my_get_expression_value(length_param_expr) *rel_shift;
      std::ostringstream ostr;
      ostr << std::setprecision(17) << at + at_shift;
      std::string expr_str=ostr.str().c_str();
      expression* new_at_expr = new_expression(expr_str.c_str(), nullptr);
      at_expr = compound_expr(new_at_expr, at + at_shift, "+", nullptr,  0,1 ); // use compound_expr to get value updated
      strcpy(at_expr->string,expr_str.c_str()); // keep the string, avoid expression with  "+ ( 0 )"
      if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " at=" << at << " at_shift=" << at_shift << " expr_str=" << expr_str << " at_expr=" << my_dump_expression(at_expr) << std::endl;
    }
    else at_expr = compound_expr(node->at_expr, at, "+", scale_expr(length_param_expr,rel_shift),  0 ,1); // use length expression and rel_shift value, this also updates the value
    place_node_at(node,to_sequ,sliced_elem,at_expr);
  }
  else
  {
    fatal_error("This is not an element ",node->name);
  }
}

static void copy_params_from_elem(command* cmd,const element* elem_inp,const std::vector<std::string> IgnoreList={})
{ // copy all active attributes from elem_inp to cmd if they exist in cmd; used to create new thin or modified thick elements
  //   ---- problematic with comments
  ElmAttr theElmAttr(elem_inp);
  if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " elem_inp->name=" << elem_inp->name << " base_type " << elem_inp->base_type->name << '\n' << theElmAttr;
  theElmAttr.TurnOff(IgnoreList);
  std::vector<std::string> active_par_to_be_used=theElmAttr.get_list_of_active_attributes();
  for(unsigned i=0;i<active_par_to_be_used.size();++i)
  {
    command_parameter* cp=return_param_recurse(active_par_to_be_used[i].c_str(),elem_inp);  //   return_param_recurse   return_param
    SetParameter_in_cmd(cmd,cp,active_par_to_be_used[i],1);
  }
}

static int set_selected_elements(el_list* the_element_list) //  modify the_element_list
{
  if (current_sequ == nullptr || current_sequ->ex_start == nullptr) // check that there is an active sequence, would otherwise crash in get_ex_range
  {
    warning("makethin selection without active sequence,", "ignored");
    return 1;
  }
  // select, flag=makethin, clear;   only resets the selection commands by  slice_select->curr=0  in mad_select.c
  // the the_element_list ist not cleared, reset it here, before evaluation of the selection for makethin
  for(int j=0; j< the_element_list->curr; ++j) // loop over the_element_list
  {
    element* el_j = the_element_list->elem[j];
    int el_j_slice_pos = name_list_pos("slice",el_j->def->par_names);
    int el_j_thick_pos = name_list_pos("thick",el_j->def->par_names); // position of thick flag      in element list
    if(el_j_slice_pos > -1) el_j->def->par->parameters[el_j_slice_pos]->double_value=1; // set all number of slices to 1
    if(el_j_thick_pos > -1) el_j->def->par->parameters[el_j_thick_pos]->double_value=0; // by default not thick
  }
  // now evaluate the selection for makethin
  node* nodes[2];  // for range check,   from nodes[0] to  nodes[1]   default is full sequence from start to end
  nodes[0] = current_sequ->ex_start;
  nodes[1] = current_sequ->ex_end;
  int slice=1; // default
  for (int i = 0; i < slice_select->curr; ++i) // loop over "select,flag=makethin" commands
  {
    name_list* nl = slice_select->commands[i]->par_names;

    const int pos_full   = name_list_pos("full", nl);
    const bool full_fl   = pos_full  > -1 && nl->inform[pos_full];  // selection with full

    const int pos_range  = name_list_pos("range", nl);
    const bool range_fl  = pos_range > -1 && nl->inform[pos_range]; // selection with range

    const int pos_slice  = name_list_pos("slice", nl);              // position of slice parameter in select command list
    const bool slice_fl  = pos_slice > -1 && nl->inform[pos_slice]; // selection with slice

    command_parameter_list* pl = slice_select->commands[i]->par;
    if (slice_fl) slice  = pl->parameters[pos_slice]->double_value; // Parameter has been read. Slice number from select command, if given overwrites slice number which may have been defined by element definition

    const int pos_thick  = name_list_pos("thick", nl);              // position of thick flag in select command list
    const bool thick_fl  = pos_slice > -1 && nl->inform[pos_thick]; // selection with slice

    if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " i=" << std::setw(2) << i  << " nl->name=" << nl->name << " full_fl=" << full_fl << " range_fl=" << range_fl << " slice_fl=" << slice_fl << " slice=" << slice << " thick_fl=" << thick_fl << '\n';
    if(full_fl) // use full sequence from start to end, the default
    {
      nodes[0] = current_sequ->ex_start;
      nodes[1] = current_sequ->ex_end;
    }
    if(range_fl)
    {
      if (current_sequ == nullptr || current_sequ->ex_start == nullptr) // check that there is an active sequence, otherwise crash in get_ex_range
      {
        warning("makethin range selection without active sequence,", "ignored");
        return 2;
      }
      if( get_ex_range(pl->parameters[pos_range]->string, current_sequ, nodes) == 0) // set start nodes[0] and end notes[1] depending on the range string
      {
        warning("empty range", "ignored");
        continue;
      }
    }
    if(range_fl) // now elements in the sequence in the range
    {
      node* c_node = nodes[0]; // for range check.  current node
      while (c_node != nullptr) // loop over nodes in range,  set slice number in elements
      {
        element* el_j = c_node->p_elem;
        name_list* nl = el_j->def->par_names;
        const int el_j_slice_pos = name_list_pos("slice",el_j->def->par_names); // position of slice parameter in element list
        const int el_j_thick_pos = name_list_pos("thick",el_j->def->par_names); // position of thick flag      in element list
        if (pass_select_el(el_j, slice_select->commands[i]) != 0) // selection on class and pattern done in pass_select. element el_j selected
        { // the element el_j passes the selection
          if(el_j_slice_pos > -1) el_j->def->par->parameters[el_j_slice_pos]->double_value=slice; // Set the element slice number to the number of slices given in the select statement.
          if(el_j_thick_pos > -1) el_j->def->par->parameters[el_j_thick_pos]->double_value=pl->parameters[pos_thick]->double_value; // Set the element thick flag to what is given in the select statement
          if(slice>1) nl->inform[el_j_slice_pos]=1;
          nl->inform[el_j_thick_pos]=1;
        } // selection
        if (c_node == nodes[1]) break; // done with last node
        c_node = c_node->next;
      } // end of while loop over nodes in range
    } // range_fl
    else // no range_fl
    {
      for(int j=0; j< the_element_list->curr; ++j) // loop over the_element_list
      {
        element* el_j = the_element_list->elem[j];
        name_list* nl = el_j->def->par_names;
        const int el_j_slice_pos = name_list_pos("slice",el_j->def->par_names);
        const int el_j_thick_pos = name_list_pos("thick",el_j->def->par_names); // position of thick flag      in element list
        if (pass_select_el(el_j, slice_select->commands[i]) != 0) // selection on class and pattern done in pass_select. element el_j selected
        { // the element el_j passes the selection
          if(el_j_slice_pos > -1) el_j->def->par->parameters[el_j_slice_pos]->double_value=slice; // Set the element slice number to the number of slices given in the select statement.
          if(el_j_thick_pos > -1) el_j->def->par->parameters[el_j_thick_pos]->double_value=pl->parameters[pos_thick]->double_value; // Set the element thick flag to what is given in the select statement
          if(slice>1 && el_j_slice_pos > -1 ) nl->inform[el_j_slice_pos]=1; // negative drift to start
          if(el_j_thick_pos> -1) nl->inform[el_j_thick_pos]=1;
        } // selection
      } // loop over the_element_list
    } // range_fl
  } // end of loop over select slice commands
  if (MaTh::Verbose) std::cout << dump_slices(the_element_list); // shows where 2 or more slices were selected
  return 0;
}

command_parameter* k0_from_angle(const command_parameter* angle_param)
{
  command_parameter* k0cmdptr=new_command_parameter("k0", k_double);
  if (angle_param->expr) k0cmdptr->expr =  clone_expression(angle_param->expr);
  k0cmdptr->double_value = angle_param->double_value;
  if ( MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << my_dump_command_parameter(k0cmdptr)  << std::endl;
  return k0cmdptr;
}

void makethin(in_cmd* incmd) // public interface to slice a sequence, called by exec_command from mad_cmd.c
{ // generates new sliced elements using my_make_element which adds to the global element_list
  double CPU_start=clock();
  name_list* nl = incmd->clone->par_names;
  command_parameter_list* pl = incmd->clone->par;

  static std::string LastSequenceSliced,LastStyle="teapot";
  static SequenceList sliced_seqlist;

  MaTh::Verbose=0;
  if(get_option("debug")) MaTh::Verbose=1;
  if(get_option("verbose")) MaTh::Verbose=2;

  const int ipos_style = name_list_pos("style", nl);
  std::string slice_style;

  if (nl->inform[ipos_style] &&  pl->parameters[ipos_style]->string )
  {
    slice_style = pl->parameters[ipos_style]->string ;
    std::cout << "makethin: style chosen : " << slice_style << '\n';
  } else slice_style = "teapot"; // Should be "hybrid" for backward compatibility

  // first check makethin parameters which influence the selection

  int iMakeConsistent;
  const int ipos_mk = name_list_pos("makeconsistent", nl);
  if( ipos_mk > -1 && nl->inform[ipos_mk])
    iMakeConsistent = pl->parameters[ipos_mk]->double_value;
  else iMakeConsistent = 0; // default is false in mad_dict.c

  const int ipos_me = name_list_pos("makeendmarkers", nl);
  if( ipos_me > -1 && nl->inform[ipos_me])
    MaTh::iMakeEndMarkers = pl->parameters[ipos_me]->double_value;
  else MaTh::iMakeEndMarkers = 0; // default is false in mad_dict.c

  const int ipos_md = name_list_pos("makedipedge", nl);
  if( ipos_md > -1 && nl->inform[ipos_md])
    MaTh::iMakeDipedge=pl->parameters[ipos_md]->double_value;
  else MaTh::iMakeDipedge = 1; // default is true in mad_dict.c

  const int ipos_mx = name_list_pos("moreexpressions", nl);
  if( ipos_mx > -1 && nl->inform[ipos_mx])
    MaTh::iMoreExpressions=pl->parameters[ipos_mx]->double_value;
  else MaTh::iMoreExpressions = 1; // default is 1 mad_dict.c

  if (MaTh::Verbose>1)
  {
    // controlled by input
    std::cout << "makethin slice_style=" << slice_style
    << "   flags: "
    << " makedipedge=" << MaTh::iMakeDipedge
    << " makeconsistent=" << iMakeConsistent
    << " makeendmarkers=" << MaTh::iMakeEndMarkers
    << " moreexpressions=" << MaTh::iMoreExpressions
    << "   hard coded flags: "
    << " kill_fringe_fl=" << kill_fringe_fl       // set to true for thick bend body slices
    << " dipedge_h1_h2_fl=" << dipedge_h1_h2_fl   // normally off
    << " Enable_all_attr_fl=" << Enable_all_attr_fl // implemented 4/2018 to allow for extra parameters on sliced sequence
    << " Verbose=" << MaTh::Verbose
    << '\n';
    // hard coded
    if(dipedge_h1_h2_fl) std::cout << "dipedge_h1_h2_fl=" << dipedge_h1_h2_fl << " is on. Higher order h1, h2 parameters will be kept. Tracking may become non-simplectic" << '\n';
  }

  if (slice_select->curr > 0)
  {
    int iret=set_selected_elements(element_list); // makethin selection -- modifies global element_list
    if (MaTh::Verbose&& iret) std::cout << "after set_selected_elements iret=" << iret << std::endl;
  }
  else  warning("makethin: no selection list,","slicing all to one thin lens.");

  if(iMakeConsistent) force_consistent_slices(element_list);

  const int ipos_seq = name_list_pos("sequence", nl);
  char* name = nullptr;
  if (nl->inform[ipos_seq] && (name = pl->parameters[ipos_seq]->string))
  {
    const int ipos2 = name_list_pos(name, sequences->list);
    if (ipos2 >= 0)
    {
      sequence* thick_sequ = sequences->sequs[ipos2];
      if(thick_sequ->ref_flag!=0)
      {
        warning("REFER in lattice must be set to CENTER, MAKETHIN:", "ignored");
      }
      else
      {
        sequence* sliced_seq = sliced_seqlist.slice_sequence(slice_style,thick_sequ,LastSequenceSliced,LastStyle); // slice the sequence
        disable_line(sliced_seq->name, line_list);
        sliced_seq->start->previous = sliced_seq->end;
        LastSequenceSliced=thick_sequ->name;
      }
    }
    else warning("unknown sequence ignored:", name);
  }
  else warning("makethin without sequence:", "ignored");
  LastStyle=slice_style;
  if (MaTh::Verbose) std::cout << "makethin: finished in " << (clock()-CPU_start)/CLOCKS_PER_SEC << " seconds" << '\n';
}

//--------  SliceDistPos
SliceDistPos::SliceDistPos(const int n,const bool teapot_fl) : delta(0.5), Delta(0),delta_str("1/2"),delta_half_str("1/4"),Delta_str("0"),Delta_half_str("0")
{ // typically called with    slice_style==std::string("teapot")    which is true for teapot  and false for simple
  // note that n = number of cuts = number of thin slices = number of thick slices -1
  // called for thick slices, positions of thin slices are calculated with simple_at_shift teapot_at_shift
  this->n=n;
  this->teapot_fl=teapot_fl;
  if(n>1)
  {
    if(teapot_fl) delta=1./(2*(1+n)); else delta=1./(2.*n);
  }
  if(n>1)
  {
    if(teapot_fl) Delta=n/(n*n-1.); else Delta=1./n;
  }

  //new same with expression
  if(n>1)
  {
    if(teapot_fl)
    {
      delta_str     ="1/"+std::to_string(2*(1+n));
      delta_half_str="1/"+std::to_string(4*(1+n));
    }
    else
    {
      delta_str     ="1/"+std::to_string(2*n);
      delta_half_str="1/"+std::to_string(4*n);
    }
  }
  if(n>1)
  {
    if(teapot_fl)
    {
      Delta_str     =std::to_string(n)+"/"+std::to_string(n*n-1);
      Delta_half_str=std::to_string(n)+"/"+std::to_string(2*(n*n-1));
    }
    else
    {
      Delta_str     ="1/"+std::to_string(n);
      Delta_half_str="1/"+std::to_string(2*n);
    }
  }
  if(MaTh::Verbose>1) Print();
}

void SliceDistPos::Print(std::ostream &StrOut) const
{
  StrOut << "SliceDistPos::Print teapot_fl=" << teapot_fl << " n=" << n << " delta=" << delta << " Delta=" << Delta
  << " delta_str=" << delta_str << " delta_half_str=" << delta_half_str << " Delta_str=" << Delta_str << " Delta_half_str=" << Delta_half_str << '\n';
}

//--------  OneElementWithSlices
OneElementWithSlices::OneElementWithSlices(const element* thick_elem,element* sliced_elem) : thick_elem(thick_elem) // OneElementWithSlices constructor
{
  theSlices.push_back(sliced_elem); // there can be several slices
}

//--------  ElementListWithSlices
ElementListWithSlices::ElementListWithSlices(unsigned int verbose) : verbose(verbose), get_thin_calls(0), get_thin_iteractions(0), ilast1(-1), ilast2(-1) // constructor, init counters
{
  if(verbose>2) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " ElementListWithSlices constructor called" << '\n';
}

ElementListWithSlices::~ElementListWithSlices() // destructor
{
  if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " ElementListWithSlices destructor called VecElemWithSlices.size()=" << VecElemWithSlices.size() << std::endl;
  for(unsigned int iel=0; iel<VecElemWithSlices.size(); ++iel) delete VecElemWithSlices[iel]; // undo the new OneElementWithSlices(thick_elem,sliced_elem);
}

void ElementListWithSlices::PrintCounter(std::ostream &StrOut) const
{
  StrOut << "ElementListWithSlices::PrintCounter " << " get_thin_calls=" << get_thin_calls << " get_thin_iteractions=" << get_thin_iteractions;
  if(VecElemWithSlices.size()>0 && get_thin_calls>0) StrOut<< " get_thin_iteractions/get_thin_calls=" << get_thin_iteractions/get_thin_calls << " ineff=" << get_thin_iteractions/((double)VecElemWithSlices.size()*get_thin_calls);
  StrOut << '\n';
}

int ElementListWithSlices::find_thick(const element* thick_elem) // find thick_element pointer in VecElemWithSlices
{
  get_thin_calls++;
  int ifound=-1;
  if(VecElemWithSlices.size()>0)
  {
    unsigned int isearched=0; // local counter of how many searched
    // look for the element in the list
    if(ilast2 > -1 && VecElemWithSlices[ilast2]->thick_elem==thick_elem)
    {
      ifound=ilast2; // same as ilast2, no search needed
    }
    else if(ilast1 > -1 && VecElemWithSlices[ilast1]->thick_elem==thick_elem)
    {
      ifound=ilast1; // same as ilast1, no search needed
    }
    else
    {
      for(unsigned int iel=0; iel<VecElemWithSlices.size(); ++iel) // search forward
        // for(int iel=VecElemWithSlices.size()-1; iel>=0; iel--)  // search backward
      {
        get_thin_iteractions++;
        isearched++; // local counter
        if( VecElemWithSlices[iel]->thick_elem == thick_elem )
        {
          ifound=iel; // thick element already known
          break;
        }
      }
    }
    if(ifound<0)
    {
      ilast2=ilast1;
      ilast1=isearched;
      if(ilast1==(int)VecElemWithSlices.size()) ilast1--; // make sure ilast remains within list
    }
    else // found
    {
      ilast2=ilast1;
      ilast1=ifound; // remember last found
    }
  }
  return ifound;
}

element* ElementListWithSlices::find_slice(const element* thick_elem,const int slice) // find address of thin slice by slice number for thick_elem
{
  const int iel=find_thick(thick_elem);
  if(iel<0)
  {
    if(verbose>1) std::cout << '\n';
    return nullptr;
  }
  // thick element found, now check if slice already defined
  int islice=slice-1;
  int nslices = (int) VecElemWithSlices[iel]->theSlices.size();

  element* result=nullptr;
  if(islice < nslices)
  {
    result=VecElemWithSlices[iel]->theSlices[islice]; // already done
  }
  else if(verbose>1) std::cout << " slice " << slice << " still to do" << '\n';
  return result;
}

element* ElementListWithSlices::find_slice(const element* thick_elem,const std::string name) // find address of thin slice by slice name for thick_elem
{
  const int iel=find_thick(thick_elem);
  if(iel<0)
  {
    return nullptr;
  }
  element* result=nullptr;
  const int nslices = (int) VecElemWithSlices[iel]->theSlices.size();
  for(unsigned int i=0; i<(unsigned int)nslices; ++i)
  {
    if( std::string(VecElemWithSlices[iel]->theSlices[i]->name)==name) // found
    {
      result=VecElemWithSlices[iel]->theSlices[i]; // can still be NULL, in case of edge elements for e1=0
      break;
    }
  }
  if(verbose>1 && result==nullptr) if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " returns nullptr for " << name << '\n';
  return result;
}

void ElementListWithSlices::put_slice(const element* thick_elem,element* sliced_elem) // add sliced_elem to the list
{
  bool found=false;
  if(thick_elem && sliced_elem)
  {
    for(unsigned int iel=0; iel<VecElemWithSlices.size(); ++iel)
    {
      if( strcmp(VecElemWithSlices[iel]->thick_elem->name,thick_elem->name) == 0 )
      {
        found=true;
        VecElemWithSlices[iel]->theSlices.push_back(sliced_elem);
        if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " put_slice found thick name=" << std::setw(MaTh::el_name_maxlen) << thick_elem->name << " slice name=" << sliced_elem->name << " in list at iel=" << iel << " #slices=" << VecElemWithSlices[iel]->theSlices.size() << '\n';
        break;
      }
    }
    if(!found)
    {
      OneElementWithSlices* aSliceList= new OneElementWithSlices(thick_elem,sliced_elem);
      VecElemWithSlices.push_back(aSliceList);
      if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " put_slice add  thick=" << std::left << std::setw(MaTh::el_name_maxlen) << thick_elem->name << std::setw(MaTh::par_name_maxlen) << " thin=" << sliced_elem->name << std::right << " to list, now VecElemWithSlices.size()=" << VecElemWithSlices.size() << '\n';
    }
  }
  else
  {
    std::ostringstream WarnStr;
    WarnStr << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " put_slice called with undefined thick_elem=" << thick_elem << " or sliced_elem=" << sliced_elem;
    warning_to_c(WarnStr);
  }
  return;
}

void ElementListWithSlices::Print(std::ostream &StrOut) const
{
  StrOut << " iel  #slices   " << std::setw(MaTh::el_type_maxlen) << "base_type" << std::setw(MaTh::par_name_maxlen) << "name" << std::setw(MaTh::par_name_maxlen) << "parent_name" << std::setw(MaTh::el_name_maxlen) << "parent->base_type" << std::setw(MaTh::par_name_maxlen) << "       slice_elem->name               slices     VecElemWithSlices.size()=" << VecElemWithSlices.size();
  StrOut << '\n';
  for(unsigned int iel=0; iel<VecElemWithSlices.size(); ++iel) // original
  {
    unsigned int nslices = (unsigned int) VecElemWithSlices[iel]->theSlices.size();
    if(verbose || nslices>1) // show only if more than 1 slice,  show all in case of verbose
    {
      const element* el_thick=VecElemWithSlices[iel]->thick_elem;
      StrOut << std::setw(4) << iel << std::setw(8) << nslices << std::setw(15) << el_thick->base_type->name << std::setw(MaTh::par_name_maxlen) << el_thick->name;
      if(el_thick && el_thick->parent)            StrOut << std::setw(MaTh::par_name_maxlen) << el_thick->parent->name; else StrOut << std::setw(MaTh::par_name_maxlen) << " ";
      if(el_thick && el_thick->parent->base_type) StrOut << std::setw(MaTh::par_name_maxlen) << el_thick->parent->base_type->name; else StrOut << std::setw(MaTh::par_name_maxlen) << " ";
      for(unsigned int i=0; i<nslices; ++i)
      {
        const element* eli=VecElemWithSlices[iel]->theSlices[i];
        if(eli) StrOut << std::setw(MaTh::par_name_maxlen) << eli->name; else StrOut << std::setw(MaTh::par_name_maxlen) << " ";
        StrOut << " address "  << std::setw(12) <<  eli;
      }
      StrOut << '\n';
    }
  }
  StrOut << '\n';
}

//--------  SeqElList
SeqElList::SeqElList(const std::string seqname,const std::string slice_style,/* sequence* thick_sequ,*/sequence* sliced_seq,node* thick_node,SequenceList* theSequenceList)
: thick_node(thick_node), work_node(thick_node), thick_elem_sliced(nullptr), knl_param(nullptr), kns_param(nullptr), theSequenceList(theSequenceList), sliced_seq(sliced_seq), seqname(seqname), slice_style(slice_style), verbose(MaTh::Verbose), eps(1.e-15), MakeDipedge(MaTh::iMakeDipedge) // SeqElList constructor, eps used to check if values are is compatible with zero
{
}

SeqElList::~SeqElList() // destructor
{
  if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << std::endl;
}

double SeqElList::simple_at_shift(const int slices,const int slice_no) const // return at relative shifts from centre of unsliced magnet
{
  const int n = slices;
  const int i = slice_no;
  return n>1 ? (2.0*i-1)/(2.0*n)-0.5 : 0.0;
}

double SeqElList::teapot_at_shift(int slices, int slice_no) const
{
  const int n = slices; // number of cuts or thin slices, gives n+1 thick slices
  const int i = slice_no;
  double shift=n>1 ? 0.5*n*(1-2*i+n)/(1.0-n*n) : 0.0; // see  http://ab-dep-abp.web.cern.ch/ab-dep-abp/LCU/LCU_meetings/2012/20120918/LCU_makethin_2012_09_18.pdf
  if ( MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " n=" << n << " i=" << i << " shift=" << shift << '\n';
  return shift;
}

double SeqElList::collim_at_shift(const int slices,const int slice_no) const
{
  const int n = slices;
  const int i = slice_no;
  return n>1 ? (i-1.0)/(n-1.0)-0.5 : 0.0;
}

double SeqElList::hybrid_at_shift(const int slices, const int slice_no) const
{
  return slices>4 ? simple_at_shift(slices, slice_no) : teapot_at_shift(slices, slice_no); // old for backwards compatibility, should be removed in future
}

double SeqElList::at_shift(const int slices,const int slice_no,const std::string local_slice_style) const// return relative shifts from centre of unsliced magnet
{
  double shift=0;
  if (!slices || !slice_no) fatal_error("makethin: invalid slicing for zero slices",local_slice_style.c_str());
  if      (local_slice_style==std::string("hybrid"))  shift= hybrid_at_shift(slices,slice_no);
  else if (local_slice_style==std::string("simple"))  shift= simple_at_shift(slices,slice_no);
  else if (local_slice_style==std::string("teapot"))  shift= teapot_at_shift(slices,slice_no);
  else if (local_slice_style==std::string("collim"))  shift= collim_at_shift(slices,slice_no);
  else fatal_error("makethin: Style chosen not known:",local_slice_style.c_str());
  if(MaTh::Verbose>1/* && local_slice_style!=slice_style */ ) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " local_slice_style=" << local_slice_style << " slice_style=" << slice_style << " shift=" << shift << '\n';
  return shift;
}

void SeqElList::Print(std::ostream& StrOut) const
{
  StrOut << "SeqElList::Print seqname=" << seqname << " theSliceList->VecElemWithSlices.size()=" << theSliceList->VecElemWithSlices.size() << " slice_style=\"" << slice_style << "\"" << '\n';

  StrOut << '\n' << "   theSliceList:" << '\n' << *theSliceList;
  if(verbose) theSliceList->PrintCounter(StrOut);

  StrOut << '\n' << "   theRbendList:" << '\n' << *theRbendList;
  if(verbose) theRbendList->PrintCounter();

  StrOut << '\n' << "theBendEdgeList:" << '\n' << *theBendEdgeList;
  if(verbose) theBendEdgeList->PrintCounter(StrOut);
}

element* SeqElList::new_marker_element(const std::string el_name, const element* el_inp)
{
  command* cmd = new_cmdptr( find_element("marker", base_type_list) );
  copy_params_from_elem(cmd,el_inp,MaTh::DoNotCopy);
  return my_El_List->my_make_element(el_name,"marker",cmd,-1); // makes new marker
}

element* SeqElList::create_wire_element(const element* thick_elem,int slice_no) // for slicing wire collimators, wire(s) with zero length
{
  element* newwire=nullptr;
  command_parameter* ilnorm_param = return_param_recurse("current",thick_elem);
  if(ilnorm_param)
  {
    command* cmd = new_cmdptr( find_element("wire", base_type_list) );
    for(unsigned i=0;i<MaTh::WireCollimatorParmList.size();++i)
    {
      command_parameter* cp=return_param_recurse(MaTh::WireCollimatorParmList[i].c_str(),thick_elem); // copy wire parameters from collimator wire
      SetParameter_in_cmd(cmd,cp,MaTh::WireCollimatorParmList[i],1);
    }
    slice_attributes_to_slice(cmd,thick_elem); // scale ilnorm with nslices
    std::string WireName=std::string(thick_elem->name)+"_wire";
    if(nslices>1) WireName=WireName+".."+std::to_string(slice_no);
    newwire= my_El_List->my_make_element(WireName,"wire",cmd,-1);
  }
  return newwire;
}

void SeqElList::kn_ks_from_thick_elem(const element* thick_elem,command_parameter* kn_pars[4],command_parameter* ks_pars[4]) const
{ // read k0-k3, k0s-k3s in thick_elem and put them in kn_pars, ks_pars;  clone to allow for changes in copy
  command_parameter* p;
  const std::vector<std::string> kn_name={"k0", "k1", "k2", "k3"};
  const std::vector<std::string> ks_name={"k0s","k1s","k2s","k3s"};
  for(unsigned int i=0;i<kn_name.size();++i) kn_pars[i] = (p=return_param_recurse(kn_name[i].c_str() ,thick_elem)) ? clone_command_parameter(p) : nullptr;
  for(unsigned int i=0;i<ks_name.size();++i) ks_pars[i] = (p=return_param_recurse(ks_name[i].c_str() ,thick_elem)) ? clone_command_parameter(p) : nullptr;
}

void SeqElList::add_ktap(command_parameter* k_param,const element* thick_elem)
{
  if(!k_param) return;
  std::string k_param_name=k_param->name;
  std::string k_name,ktap_name;
  for(unsigned int i=1;i<3;++i)
  {
    k_name="k"+std::to_string(i); // k1, k2
    if(k_param_name=="ksl") k_name=k_name+"s"; // or  k1s, k2s
    ktap_name=k_name+"tap"; // corresponding tap attribute
    add_ktap_i(i,k_param,k_name,ktap_name,thick_elem);
  }
}

void SeqElList::add_ktap_i(const int i,command_parameter* k_param,const std::string k_name,const std::string ktap_name,const element* thick_elem)
{
  if( command_parameter *p     = return_param_recurse(k_name.c_str(), thick_elem) )  // has k1n, or k1s, k2n, k2s
    if( command_parameter *p_tap = return_param_recurse(ktap_name.c_str(), thick_elem) )// has corresponding tap
    {
      if(p->expr) k_param->expr_list->list[i] = compound_expr(p->expr,                          p->double_value, "+", p_tap->expr, p_tap->double_value, 0);
      else        k_param->expr_list->list[i] = compound_expr(expr_from_value(p->double_value), p->double_value, "+", p_tap->expr, p_tap->double_value, 0);
    }
}

command_parameter* SeqElList::make_k_list(const std::string parnam,command_parameter* k_pars[4]) const
{
  command_parameter* k_param=nullptr;
  // from k values 0-3 to  expr lists with these values
  if ( k_pars[0] || k_pars[1] || k_pars[2] || k_pars[3] )
  {
    k_param = new_command_parameter(parnam.c_str(), k_double_array);
    k_param->expr_list = new_expr_list(10);
    k_param->double_array = new_double_array(10);

    for(int i=0; i<4; ++i) // set k_param with expr_list
    { // initialize the parameters with nullptr for expression and 0 for the value
      k_param->expr_list->list[i] = nullptr;
      k_param->double_array->a[i] = 0;
      if (k_pars[i]) // copy existing k's
      {
        if (k_pars[i]->expr) k_param->expr_list->list[i] = clone_expression(k_pars[i]->expr);
        k_param->double_array->a[i] = k_pars[i]->double_value;
      }
      k_param->expr_list->curr++; k_param->double_array->curr++; // update the number of k's in our arrays
    }
  }
  return k_param;
}

element* SeqElList::create_bend_dipedge_element(element* thick_elem,const bool Entry) // using Dipedge  http://mad.web.cern.ch/mad/madx.old/Introduction/dipedge.html
{
  // see also twiss.f90    SUBROUTINE tmbend,   and  SUBROUTINE tmfrng  for fringe fields, and   SUBROUTINE tmdpdg  for dipedge
  // makes dipedge element for start or end of the dipole
  // example
  // from original thick
  // mb1: sbend,l:=lmb ,angle:=amb ,k1:=kmb ,e1:=ee1 ,e2:=ee2 ;
  // to
  // mb1_l: dipedge, h:= amb/lmb ; e1:=ee1 ;  !--------- new dipedge at entry
  // mb1: sbend,l:=lmb ,angle:=amb ,k1:=kmb ;      ! modified middle, e1, e2  removed
  // mb1_r: dipedge, h:= amb/lmb ; e1:=ee2 ;  !--------- new dipedge at exit
  //
  // request from Laurent Deniau and Andrea Latina in 10/2014   also move any h1, h2  parameters as h parameter to entry, exit dipedge
  //

  static std::vector<std::string> CheckBendParams = {
    "polarity", "tilt", "hgap", "mech_sep", "v_pos", "magnet", "model", "method", "exact", "nst" };

  element* dipedge=NULL;
  std::string thick_elem_name=thick_elem->name;
  if(thick_elem_name[0]==MaTh::ExtraChar)
  {
    thick_elem_name=thick_elem_name.substr(1); // without the ExtraChar
  }

  if (thick_elem)
  {
    std::string dipedge_name=thick_elem_name;
    if(Entry) dipedge_name+="_den"; else dipedge_name+="_dex"; // dipedge entry or exit

    std::string dipedge_cmd_name="dipedge";
    if(Entry) dipedge_cmd_name+="_l_"; else dipedge_cmd_name+="_r_";
    dipedge_cmd_name+="cmd";

    const double eps=1.e-15;

    expression* l_par_expr=my_get_param_expression(thick_elem, "l"); // with this l_par_expr should not be NULL
    expression* angle_par_expr = my_get_param_expression(thick_elem,"angle");
    command_parameter* hparam=new_command_parameter("h",k_double);
    hparam->expr=compound_expr(angle_par_expr,0.,"/",l_par_expr,0,1); // this also updates the value

    command* dipedge_cmd = new_cmdptr( find_element("dipedge", base_type_list) );

    SetParameter_in_cmd(dipedge_cmd,hparam,"h",1);
    if(dipedge_h1_h2_fl)
    {
      if(Entry) SetParameter_in_cmd(dipedge_cmd, return_param_recurse("h1",thick_elem), "h1",1); // at entry, copy h1 from thick bend as dipedge h1
      else      SetParameter_in_cmd(dipedge_cmd, return_param_recurse("h2",thick_elem), "h2",1);  // at  exit, copy h2 from thick bend as dipedge h2
    }
    if(Entry)   SetParameter_in_cmd(dipedge_cmd, return_param_recurse("e1",thick_elem), "e1",1); // at entry, copy e1 from thick sbend as dipedge e1
    else        SetParameter_in_cmd(dipedge_cmd, return_param_recurse("e2",thick_elem), "e1",1); // at  exit, copy e2 from thick sbend as dipedge e1

    if(Entry)
    {
      SetParameter_in_cmd(dipedge_cmd, return_param_recurse("fint",thick_elem), "fint",1);
    }
    else // Exit
    {
      const command_parameter *fintxparam = return_param_recurse("fintx",(const element*)thick_elem); // use my const version of return_param to check the presence of fintx in thick_elem
      const command_parameter *fintparam  = return_param_recurse("fint", (const element*)thick_elem);
      if(fintxparam!=nullptr) SetParameter_in_cmd(dipedge_cmd, fintxparam, "fint",1); // use fintx if given, no check on value, allows for fintx=0 to turn off exit
      else                    SetParameter_in_cmd(dipedge_cmd, fintparam,  "fint",1); // use fint
      if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << my_dump_command_parameter(fintxparam) << '\n';
    }

    if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << my_dump_command(dipedge_cmd) << '\n';

    bool found=false;

    for(unsigned int i=0; i<CheckBendParams.size(); ++i) // copy other nontrivial parameters given in CheckParams from thick bend  -- only gets here with nontrivial e1 or e2     -- otherwise already returned NULL before
    {
      const std::string parnam = CheckBendParams[i];
      command_parameter* this_param = return_param_recurse(parnam.c_str(), thick_elem);
      double value      =my_get_int_or_double_value(thick_elem           ,parnam,found);
      double default_val=my_get_int_or_double_value(thick_elem->base_type,parnam,found);
      if(this_param || (found && fabs(value-default_val)>eps) ) // expression defined or non-trivial value
      {
        SetParameter_in_cmd(dipedge_cmd,this_param,parnam,1);
        if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << thick_elem->name << "        use parameter " << std::setw(12) << parnam << " for dipedge this_param=" << this_param << '\n';
      }
    }

    dipedge=my_El_List->my_make_element(dipedge_name, "dipedge", dipedge_cmd,-1); // make the dipedge element using the command dipedge_cmd, -1 means avoid warnings,  1 means delete and warn, 2 means warn and ignore if already present  see  add_to_el_list  mad_elem.c

  }
  return dipedge;
}

element* SeqElList::sbend_from_rbend(element* rbend_el)
{
  bool has_angle;
  double angle=my_get_int_or_double_value(rbend_el,"angle",has_angle);
  if(has_angle && fabs(angle)<eps) has_angle=false;
  if(!has_angle) return rbend_el; // angle needed to translate length


  element* sbend_el_parent;
  if (rbend_el == rbend_el->parent) return nullptr; // no further parent to consider
  else
  {
    element* el_found=theRbendList->find_slice(rbend_el->parent,MaTh::ExtraChar+std::string(rbend_el->parent->name)); // check if parent already translated
    if(!el_found)
    {
      sbend_el_parent = sbend_from_rbend(rbend_el->parent); // recursively translate parent
      if( verbose>1 ) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " after  recursive call sbend_from_rbend sbend_el_parent->name=" << sbend_el_parent->name << std::endl;
    }
  }

  const std::string rbend_name=rbend_el->name;
  const std::string sbend_name=MaTh::ExtraChar+rbend_name; // add ExtraChar to be able to distinguish internally, avoids test-track-10  bend with zero length: bw2.ql12.r1, for slices use  thick_elem_name without the extra MaTh::ExtraChar

  command* sbend_cmd = new_cmdptr( find_element("sbend", base_type_list) );
  copy_params_from_elem(sbend_cmd,rbend_el,MaTh::DoNotCopy2);

  // --   at this point, the rbend length should already be in sbend_cmd, but needs to be modified when rbarc is on (default)
  if(rbarc_fl())
  {
    expression* l_sbend_expr=curved_from_straight_length(rbend_el); // use this modified length expression in sbend_cmd
    int il=name_list_pos("l",sbend_cmd->par_names); // parameter 0
    if(il > -1) sbend_cmd->par->parameters[il]->expr=l_sbend_expr;
    if( verbose>1 ) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " after increase of rbend length now l_sbend_expr : " << my_dump_expression(l_sbend_expr) << '\n';
  }

  element* sbend_el=my_El_List->my_make_element(sbend_name, "sbend", sbend_cmd,-1);

  if(rbend_el->parent)
  { // inherit hierarchy
    sbend_el->parent=rbend_el->parent;
  }

  add_half_angle_to(rbend_el,sbend_el,"e1");
  add_half_angle_to(rbend_el,sbend_el,"e2");

  theRbendList->put_slice(rbend_el,sbend_el); // keep address of translated rbend_el

  return sbend_el;
} // sbend_from_rbend

element* SeqElList::create_thick_slice(const element* thick_elem,const int slice_type) // create entry/body/exit slice elements
{
  const int n=nslices-1; // in case of thick slices one less
  const char* eltype=thick_elem->base_type->name;
  std::string slice_name;
  const bool entry_fl = slice_type == 0;
  const bool  exit_fl = slice_type == 2;
  const bool IsBend = strcmp(work_node->base_name, "sbend") == 0;
  std::string thick_elem_name=thick_elem->name;
  if(thick_elem_name[0]==MaTh::ExtraChar)
  {
    thick_elem_name=thick_elem_name.substr(1); // without the ExtraChar
  }
  if      (entry_fl)  slice_name=thick_elem_name+"_en"; // entry
  else if (exit_fl)   slice_name=thick_elem_name+"_ex"; // exit
  else                slice_name=thick_elem_name+"_bo"; // body  slice

  element* sliced_elem;
  if( (sliced_elem = theSliceList->find_slice(thick_elem,slice_name)))
  {
    if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " slice_name already exists, use it" << '\n';
    return sliced_elem;
  }

  SliceDistPos SP(n, slice_style==std::string("teapot") );
  if(verbose>1)
  {
    std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << eltype << " create " << slice_name << " based on " << thick_elem->name
    << " slice_type=" << slice_type << " n=" << n
    << " entry_fl=" << entry_fl
    << " exit_fl="  << exit_fl;
  }

  expression*       l_par_expr=my_get_param_expression(thick_elem, "l");     // get length expression
  expression* angle_par_expr = my_get_param_expression(thick_elem, "angle"); // get angle expressions - relevant for bends

  if(l_par_expr==nullptr) // compound_expr with scaling will fail   -- should never happen
  {
    std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " *** error *** l_par_expr=" << l_par_expr << '\n';
    exit(1);
  }

  double LengthFraction;
  if(entry_fl || exit_fl) LengthFraction=SP.delta; // start / end slice
  else                    LengthFraction=SP.Delta; // the middle or body piece

  l_par_expr                        = compound_expr(l_par_expr,    0, "*", nullptr, LengthFraction,1); // multiply length parameter expression with LengthFraction
  if(angle_par_expr) angle_par_expr = compound_expr(angle_par_expr,0, "*", nullptr, LengthFraction,1); // multiply angle  parameter expression with LengthFraction, only relevant for bends

  command* cmd = new_cmdptr( thick_elem );
  copy_params_from_elem(cmd,thick_elem,MaTh::DoNotCopy);

  if(verbose>1)
  {
    std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " thick_elem " << thick_elem->name << " " << SP;
    std::cout  << " LengthFraction=" << LengthFraction << " scaled l_par_expr " << my_dump_expression(l_par_expr);
    if(angle_par_expr) std::cout << " scaled angle_par_expr " << my_dump_expression(angle_par_expr);
    std::cout << '\n';
  }

  // now use the scaled length and if relevant angle parameters to set up the new sliced_elem via cmd
  const int length_i = name_list_pos("l", cmd->par_names);
  if(length_i > -1)
  {
    cmd->par->parameters[length_i]->expr=l_par_expr; // use the length expression in cmd to set up sliced_elem
  }
  else
  {
    std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " *** error *** thick_elem " <<  thick_elem->name << " has no length parameter : " << my_dump_element(thick_elem);
    return nullptr;
  }

  const int angle_i = name_list_pos("angle",cmd->par_names);
  if(verbose>1) std::cout << '\n' << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " angle_i=" << angle_i << '\n';
  if(angle_i > -1)
  {
    cmd->par->parameters[angle_i]->expr=angle_par_expr; // use the scaled angle_par_expr in cmd to set up sliced_elem
  }

  sliced_elem = my_El_List->my_make_element(slice_name, eltype,cmd,-1); // make new element (without parent) using the command cmd, -1 means avoid warnings

  if(length_i > -1) ParameterTurnOn("l",     sliced_elem);
  if(angle_i  > -1) ParameterTurnOn("angle", sliced_elem);
  ParameterRemove("slice", sliced_elem); // slicing done, no reason to leave the slice parameter
  ParameterTurnOn("thick", sliced_elem); //-- so that thick=true is written  to signal this for thick tracking

  if(IsBend)
  {
    if(entry_fl)
    { // bend entry,  remove/turn off exit parameters
      ParameterRemove("e2"   ,sliced_elem);
      ParameterRemove("h2"   ,sliced_elem);
      ParameterTurnOn("fint" ,sliced_elem);  // use fint
      ParameterTurnOn("fintx",sliced_elem);
      SetParameterValue("fintx",sliced_elem,0); // leave fintx, with value 0, otherwise taking fint
      SetParameterValue("kill_exi_fringe",sliced_elem,true,k_logical);
      ParameterTurnOn("kill_exi_fringe"  ,sliced_elem); // turn writing on
    }
    else if(exit_fl)
    { // bend exit, remove entry parameters
      if (verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << '\n';
      ParameterRemove("e1",sliced_elem);
      ParameterRemove("h1",sliced_elem);
      SetParameterValue("kill_ent_fringe",sliced_elem,true,k_logical);
      ParameterTurnOn("kill_ent_fringe"  ,sliced_elem); // turn writing on
      const command* el_def=sliced_elem->def;
      name_list* nl=el_def->par_names;
      int i_fint = name_list_pos("fint", nl);
      const bool fint_on =  (i_fint > -1) && el_def->par_names->inform[i_fint];
      int i_fintx = name_list_pos("fintx", nl);
      const bool fintx_on =  (i_fintx > -1) && el_def->par_names->inform[i_fintx];
      if(fintx_on) ParameterRemove("fint", sliced_elem); // fintx is on and will be used, just remove any fint on the exit
      else if(fint_on)
      { // no fintx, use fint as fintx for exit
        ParameterTurnOn("fintx",sliced_elem);
        if(i_fintx) // should be there, just inform off
        {
          bool found=false;
          double fint_value=my_get_int_or_double_value(sliced_elem,"fint",found);
          SetParameterValue("fintx",sliced_elem,fint_value);
          if (verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " no fintx, use fint value " << fint_value << " as fintx for exit" << '\n';
        }
        ParameterRemove("fint",sliced_elem); // remove fint on exit
      }
    }
    else Remove_All_Fringe_Field_Parameters(sliced_elem); // thick magnet body, remove fringe fields
  }
  theSliceList->put_slice(thick_elem,sliced_elem); //-- store what is done in theSliceList
  return sliced_elem;
} // create_thick_slice

element* SeqElList::create_thin_slices(const element* thick_elem, int slice_no) // create thin magnet slices using multipole, recursively for parents
{
  element *sliced_elem_parent;
  if (thick_elem == thick_elem->parent) return nullptr; // no further parent to consider
  else sliced_elem_parent = create_thin_slices(thick_elem->parent,slice_no); // recursively slice parent       --- element specific
  std::string parent_or_base="multipole";
  if(sliced_elem_parent) parent_or_base=sliced_elem_parent->name;
  if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << std::left << std::setw(MaTh::el_name_maxlen) << thick_elem->name << " parent_or_base= " << std::setw(MaTh::el_type_maxlen) << parent_or_base << " thick_elem base_type=" << thick_elem->base_type->name << " slice_no=" << slice_no << " nslices=" << nslices << '\n';
  element *sliced_elem = theSliceList->find_slice(thick_elem,slice_no);
  if( sliced_elem )
  {
    return sliced_elem; // already done
  }
    std::string parent_name;
    if(sliced_elem_parent) parent_name=sliced_elem_parent->name;

  // work on the parameters, specific for magnets
  command_parameter* length_param = return_param_recurse("l",thick_elem);
  const command_parameter* angle_param  = return_param_recurse("angle",thick_elem);

  if( std::string(thick_elem->base_type->name)=="rbend")
  {
    if(rbarc_fl())
    {
      if(length_param && angle_param)
      {
        length_param->expr=curved_from_straight_length(thick_elem);
        length_param->double_value=my_get_expression_value(length_param->expr);
      }
    }
  }

  command_parameter *kn_pars[4],*ks_pars[4];
  kn_ks_from_thick_elem(thick_elem,kn_pars,ks_pars); // get kn_pars[0-3] and ks_pars[0-3] from thick_elem, as first step to construct kn_l ks_l lists

  bool mult_with_length = true;
  command_parameter* multipole_angle_param = nullptr;

  if(angle_param) // use angle to generate k0, if k0 non existing
  {
    double k0val=0,angleval=0,lenval=0;
    if(kn_pars[0]) // angle_param and k0 defined
    {
      if(length_param)
      {
        if(length_param->expr) lenval=my_get_expression_value(length_param->expr);
        else lenval=length_param->double_value;
      }
      if(angle_param->expr) angleval=my_get_expression_value(angle_param->expr);
      else angleval=angle_param->double_value;
      if(kn_pars[0]->expr) k0val=my_get_expression_value(kn_pars[0]->expr);
      else k0val=kn_pars[0]->double_value;
      if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << std::left << std::setw(MaTh::el_name_maxlen) << thick_elem->name << " angleval=" << angleval << " k0val=" << k0val << " lenval=" << lenval << " k0val*lenval=" << k0val*lenval << std::endl;
      if(fabs(k0val*lenval-angleval)>eps && fabs(k0val) > eps)
      { // both angle and k0 given with different information
        if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << std::left << std::setw(MaTh::el_name_maxlen) << thick_elem->name << " angle and k0 both given and not equal, write also angle to multipole" << std::endl;
        multipole_angle_param=new_command_parameter("angle", k_double);
        if (angle_param->expr) multipole_angle_param->expr = clone_expression(angle_param->expr);
        multipole_angle_param->double_value = angle_param->double_value;
      }
      if(fabs(k0val) < eps)
      { // This is in case k0 is defined but 0
        kn_pars[0]=k0_from_angle(angle_param); // k0 generated from angle
        mult_with_length = false; // angle is already  k0 * length
        if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << std::left << std::setw(MaTh::el_name_maxlen) << thick_elem->name << " k0 is defined but 0  kn_pars[0]=" << kn_pars[0] << std::endl;
      }
    }
    else // k0 not defined, angle given, generate k0
    {
      kn_pars[0]=k0_from_angle(angle_param); // k0 generated from angle
      mult_with_length = false; // angle is already  k0 * length
    }
    if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " has angle_param,  kn_pars[0] " << my_dump_command_parameter(kn_pars[0]) << std::endl;
  }

  if(thick_elem != thick_elem_sliced)
  {
    knl_param=make_k_list("knl",kn_pars); // make new knl 0,1,2,3 expr_list based on kn_pars from thick
    kns_param=make_k_list("ksl",ks_pars); // make new ksl 0,1,2,3 expr_list based on ks_pars from thick

    add_ktap(knl_param,thick_elem);
    add_ktap(kns_param,thick_elem);

    // multiply the k by length and divide by slice
    knl_param = scale_and_slice(knl_param,length_param,nslices,mult_with_length);
    kns_param = scale_and_slice(kns_param,length_param,nslices,mult_with_length);
  }

  if(multipole_angle_param && nslices>1) // case of angle different from k0l
  {
    if (angle_param->expr) multipole_angle_param->expr  =  compound_expr(angle_param->expr,0.,"/",nullptr,nslices,1);
    else multipole_angle_param->double_value /= nslices;
  }

  command* cmd = new_cmdptr( find_element("multipole", base_type_list) );
  copy_params_from_elem(cmd,thick_elem,str_v_join(MaTh::DoNotCopy,{"angle","kill_ent_fringe","kill_exi_fringe"}));
  SetParameter_in_cmd(cmd,knl_param,"knl",1);
  SetParameter_in_cmd(cmd,kns_param,"ksl",1);
  if(multipole_angle_param) SetParameter_in_cmd(cmd,multipole_angle_param,"angle",1);
  set_lrad(cmd,length_param,nslices); // keep l  as lrad
  finish_make_sliced_elem(sliced_elem, thick_elem, cmd, parent_or_base, slice_no);
  ParameterRemove("l",sliced_elem);
  thick_elem_sliced=thick_elem;
  return sliced_elem;
}

element* SeqElList::create_thin_solenoid(const element* thick_elem, int slice_no) // create the sliced element, recursively for parents
{
  element *sliced_elem_parent;
  if (thick_elem == thick_elem->parent) return nullptr; // no further parent to consider
  else sliced_elem_parent = create_thin_solenoid(thick_elem->parent,slice_no); // recursively slice parent       --- element specific
  std::string parent_or_base= thick_elem->base_type->name;
  if(sliced_elem_parent) parent_or_base=sliced_elem_parent->name;
  if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << std::left << std::setw(MaTh::el_name_maxlen) << thick_elem->name << " parent_or_base= " << std::setw(MaTh::el_type_maxlen) << parent_or_base << " slice_no=" << slice_no << " nslices=" << nslices << '\n';
  element *sliced_elem;
  if( (sliced_elem = theSliceList->find_slice(thick_elem,slice_no)) ) return sliced_elem; // already done

  // work on the parameters, specific for solenoid
  const command_parameter* length_param  = return_param_recurse("l",thick_elem);
  const command_parameter* kns_param     = return_param_recurse("ks",thick_elem);
  command_parameter* ksi_par=par_scaled(kns_param,length_param,"ksi",nslices);
  command* cmd = new_cmdptr( thick_elem );
  copy_params_from_elem(cmd,thick_elem,str_v_join(MaTh::DoNotCopy,{"kill_ent_fringe","kill_exi_fringe"}));
  SetParameter_in_cmd(cmd,ksi_par,"ksi",1);

  set_lrad(cmd,length_param,nslices); // keep l  as lrad
  if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << my_dump_command(cmd) << std::endl;
  finish_make_sliced_elem(sliced_elem, thick_elem, cmd, parent_or_base, slice_no);
  ParameterRemove("l",sliced_elem);
  return sliced_elem;
}

element* SeqElList::create_thin_elseparator(const element* thick_elem, int slice_no) // create thin elseparator element, similar to create_thin_solenoid
{
  element *sliced_elem_parent;
  if (thick_elem == thick_elem->parent) return nullptr; // no further parent to consider
  else sliced_elem_parent = create_thin_elseparator(thick_elem->parent,slice_no); // recursively slice parent       --- element specific
  std::string parent_or_base= thick_elem->base_type->name;
  if(sliced_elem_parent) parent_or_base=sliced_elem_parent->name;
  if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << std::left << std::setw(MaTh::el_name_maxlen) << thick_elem->name << " parent_or_base= " << std::setw(MaTh::el_type_maxlen) << parent_or_base << " slice_no=" << slice_no << " nslices=" << nslices << '\n';
  element *sliced_elem;
  if( (sliced_elem = theSliceList->find_slice(thick_elem,slice_no)) ) return sliced_elem; // already done

  // get parameters from the thick elseparator element
  const command_parameter* length_param  = return_param_recurse("l",thick_elem);
  const command_parameter* ex_param      = return_param_recurse("ex",thick_elem);
  const command_parameter* ey_param      = return_param_recurse("ey",thick_elem);
  command_parameter* ex_l_par=par_scaled(ex_param,length_param,"ex_l",nslices);
  command_parameter* ey_l_par=par_scaled(ey_param,length_param,"ey_l",nslices);
  command* cmd = new_cmdptr( thick_elem );
  copy_params_from_elem(cmd,thick_elem,str_v_join(MaTh::DoNotCopy,{"ex","ey"}));
  SetParameter_in_cmd(cmd,ex_l_par,"ex_l",1);
  SetParameter_in_cmd(cmd,ey_l_par,"ey_l",1);
  set_lrad(cmd,length_param,nslices); // keep l  as lrad

  if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << my_dump_command(cmd) << std::endl;
  finish_make_sliced_elem(sliced_elem, thick_elem, cmd, parent_or_base, slice_no);

  ParameterRemove("l",sliced_elem);
  return sliced_elem;
}

void SeqElList::slice_node_translate() // slice/translate and add slices to sliced sequence  sliced_seq
{
  bool UseDipedges=true; // Normally true.   For tests set false, to see the result without dipedges, dipedges will then be created but not written


  std::string base_name=thick_node->base_name;
  work_node=thick_node;
  const element* original_rbend=nullptr;
  if( base_name=="rbend")
  {
    original_rbend=thick_node->p_elem;   // original rbend, used only in search if already translated
    work_node=clone_node(thick_node, 0); // work on node copy that can be changed to sbend leaving the original
  }
  element* thick_elem=work_node->p_elem;     // work directly on this element  --  to do translation only once, element maybe used several times in sequence

  const bool ThickSLice=thick_fl(thick_elem);
  const bool IsQuad     = base_name=="quadrupole";
  const bool IsSolenoid = base_name=="solenoid";
  const bool IsBend     = base_name=="rbend" || base_name=="sbend";

  std::string thick_elem_name=thick_elem->name;
  if(thick_elem_name[0]==MaTh::ExtraChar)
  {
    thick_elem_name=thick_elem_name.substr(1); // without the ExtraChar
  }


  bool already_translated_rbend=false;
  if( base_name=="rbend" && nslices>0 ) // rbend to sbend
  {
    if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " search for " << thick_elem->name << std::endl;
    if( element* sbend_el = theRbendList->find_slice(original_rbend,thick_elem->name) )  //  see if already translated   search for thick_elem->name with ExtraChar
    {
      already_translated_rbend=true;
      thick_elem=sbend_el;           // continue with already translated sbend
    }
    else
    { // do rbend to sbend translation
      thick_elem=sbend_from_rbend(thick_elem); // translate any rbend (with both length and angle) to sbend
      strcpy(work_node->base_name,"sbend"); // set also the node base_name to sbend
    }
    work_node->p_elem=thick_elem;  // update node
    if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << "  after translate update work_node " << my_dump_node(work_node) << '\n';
  } // done with any rbend, for the rest all bends will be sbend

  if(ThickSLice && nslices>1 && !IsQuad && !IsBend && !IsSolenoid)
  {
    if(nslices>1)
    {
      std::ostringstream WarnStr;
      WarnStr << thick_elem->name << " nslices=" << nslices << " thick slicing with nslices>1 for " << work_node->base_name << " not defined. Will use nslices=1";
      warning_to_c(WarnStr);
    }
    nslices=1;
  }

  if(nslices<1 || (nslices==1 && ThickSLice && !IsBend) )
  {
    if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " node " << work_node->name << " " << thick_elem->name << " ThickSLice=" << ThickSLice << " nslices=" << nslices << " place node " << work_node->name << " without slicing" << '\n';
    add_node_at_end_of_sequence(work_node,sliced_seq); // straight copy
    return;
  }

  //--- done with initial checks, prepare for actual slicing

  if(verbose>1)
  {
    std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << thick_elem->name << " " << work_node->base_name << " slice_style=\"" << slice_style << "\""
    << " nslices=" << nslices << " IsQuad=" << IsQuad << " IsSolenoid=" << IsSolenoid << " IsBend=" << IsBend << " ThickSLice=" << ThickSLice << " at_value=" << std::setw(10) << work_node->at_value << " from_name=";
    if(work_node->from_name) std::cout << work_node->from_name; else std::cout << "NULL ";
  }

  element *EntryDipedge=nullptr, *ExitDipedge=nullptr;
  element *en = nullptr, *bo = nullptr, *ex = nullptr; // pointers to thick  entry, body, exit

  if(IsBend && MakeDipedge)
  {
    std::string thick_elem_name=thick_elem->name;
    if(thick_elem_name[0]==MaTh::ExtraChar)
    {
      thick_elem_name=thick_elem_name.substr(1); // without the ExtraChar
    }
    // find any existing EntryDipedge, sbend_el, ExitDipedge    and use them
    EntryDipedge=theBendEdgeList->find_slice(thick_elem,thick_elem_name+"_den"); // dipedge entry, NULL if not yet known or e1=0
    ExitDipedge =theBendEdgeList->find_slice(thick_elem,thick_elem_name+"_dex"); // dipedge exit,  BULL if not yet known or e2=0
    if(verbose>1)
    {
      if(verbose>1)    std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__  << "              " << std::setw(MaTh::par_name_maxlen) <<   thick_elem->name << " " << work_node->base_name << '\n';
      if(EntryDipedge) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__  << " EntryDipedge=" << std::setw(MaTh::par_name_maxlen) << EntryDipedge->name << " already exists " << EntryDipedge << '\n';
      if(ExitDipedge)  std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__  << "  ExitDipedge=" << std::setw(MaTh::par_name_maxlen) <<  ExitDipedge->name << " already exists " << EntryDipedge << '\n';
    }

    if(EntryDipedge==nullptr) // create new EntryDipedge for this bend
    { // first look if e1 or h1 are there
      const command_parameter   *e1param = return_param_recurse("e1"  ,(const element*) thick_elem);
      if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " e1param=" << e1param << " cmd_par_val(e1param)=" << cmd_par_val(e1param) << '\n';
      if (command_par_value("kill_ent_fringe",thick_elem->def) == false)
      {
        EntryDipedge=create_bend_dipedge_element(thick_elem,true); // make new StartEdge element and remove e1 from thick_elem
        theBendEdgeList->put_slice(thick_elem,EntryDipedge);       // to remember this has been translated
        if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__  << " thick node " << work_node->name << " " << thick_elem->name << " now has EntryDipedge=" << EntryDipedge << '\n';
      }
    }

    if(ExitDipedge==nullptr) // create new ExitDipedge for this bend
    {
      if (command_par_value("kill_exi_fringe",thick_elem->def) == false)
      {
        ExitDipedge=create_bend_dipedge_element(thick_elem,false); // make new ExitDipedge element and remove e2 from thick_elem
        theBendEdgeList->put_slice(thick_elem,ExitDipedge);        // to remember this has been translated
      }
    } // new ExitDipedge
    if(!already_translated_rbend) Remove_All_Fringe_Field_Parameters(thick_elem); // remove what is taken care of by dipedges
  } // IsBend && MakeDipedge

  command_parameter* length_param = return_param_recurse("l",thick_elem); // get length, value or expression
  std::string local_slice_style=slice_style; // work here with a copy of slice_style that can be modified for collim

  if(ThickSLice) // create entry, body, exit pieces,  for bends, quadrupoles, solenoids  --- if not yet existing
  {
    if(nslices==1) // single thick
    {
      en = thick_elem; // full slice as entry, no body/exit
      if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " ThickSLice, nslices=" << nslices << " create thick slices _en, _bo, _ex thick_elem->name=" << thick_elem->name << " here single slice just entry, not body, exit" << '\n';
    }
    else // nslices>1
    {
      if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " ThickSLice, nslices=" << nslices << " create thick slices _en, _bo, _ex thick_elem->name=" << thick_elem->name << '\n';
      en =               create_thick_slice(thick_elem,0); // entry slice
      if(nslices>2) bo = create_thick_slice(thick_elem,1); // body slices,   last parameter = 1     since there will be only one type of body
      ex =               create_thick_slice(thick_elem,2); // exit slice
    }
  }

  expression* at_expr = work_node->at_expr;
  double at = work_node->at_value;
  if(at_expr==nullptr) at_expr=expr_from_value_2(at); // make a new expression from the value

  double length = work_node->length; // direct curved thick_elem->length

  expression* l_expr = nullptr;
  if (length_param) l_expr  = length_param->expr;

  int middle=-1;
  if (nslices>1) middle = nslices/2+1; // used to determine after which slide to place the central marker in case of thin slicing

  if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " length_param=" << length_param << " l_expr=" << l_expr << " work_node->p_elem=" << work_node->p_elem << " thick_elem=" << thick_elem << '\n';
  if(EntryDipedge && UseDipedges) place_thin_slice(work_node,sliced_seq,EntryDipedge,-0.5); // subtract half of the length to be at start

  if(MaTh::iMakeEndMarkers) place_start_or_end_marker(true); // start

  expression* thin_at_expr=nullptr;

  for (int i=1; i<=nslices; ++i) // loop to place the nslices in the sequence
  {
    element *slice_i=nullptr;

    if(ThickSLice) // fill space between slices
    {
      if(i==1)           place_thick_slice(thick_elem,en,i); // place entry  slice
      else if(i<nslices) place_thick_slice(thick_elem,bo,i); // place body/middle slice
      else               place_thick_slice(thick_elem,ex,i); // place exit slice
      // place exit body after loop
    }
    else
    {
      if (strstr(work_node->base_name,"elseparator")) slice_i = create_thin_elseparator(thick_elem,i);
      else // magnet which can be sliced
      {
        if (IsSolenoid && !ThickSLice) slice_i = create_thin_solenoid(thick_elem,i); // create the solenoid slices
        else slice_i = create_thin_slices(thick_elem,i); // create the magnet slices
      }
      if (fabs(at_shift(nslices,i,local_slice_style))>0.0)
      {
        if( MaTh::iMoreExpressions<1 ) thin_at_expr = compound_expr(at_expr,work_node->at_value,"+",  nullptr, length*at_shift(nslices,i,local_slice_style) ,1);  // use length and shift values, no expressions
        else thin_at_expr = compound_expr(at_expr,work_node->at_value,"+",scale_expr(l_expr,at_shift(nslices,i,local_slice_style)),length*at_shift(nslices,i,local_slice_style),1); // use length expressions
      }
      else
      {
        if (at_expr) thin_at_expr = clone_expression(at_expr);
      }
    }

    if(verbose>1)
    {
      std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " work_node->base_name=" << " el->name=" << std::left << std::setw(MaTh::el_type_maxlen) << work_node->base_name << " i=" << i << " middle=" << middle << " nslices=" << nslices;
      if(slice_i) std::cout << " slice_i->name=" << slice_i->name << " at_expr " << my_dump_expression(at_expr) << " l_expr " << my_dump_expression(l_expr);
      std::cout << '\n';
    }

    if (i==middle && !ThickSLice) // create and place new marker with name of thick_elem  in the middle = poition of thick_elem
    {
      element* middle_marker=new_marker_element(thick_elem_name,thick_elem);
      place_node_at(work_node, sliced_seq, middle_marker,at_expr);
      if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " work_node " << work_node->name << " nslices=" << nslices << " i=" << i << " middle=" << middle << " place middle marker at=" << work_node->at_value << " at_expr=" << at_expr
        << " thin_at_value=" << my_get_expression_value(thin_at_expr) << '\n';
    }
    if(slice_i && !ThickSLice) place_node_at(work_node, sliced_seq, slice_i,thin_at_expr);
  }

  if(MaTh::iMakeEndMarkers) place_start_or_end_marker(false); // end

  if(ExitDipedge && UseDipedges) place_thin_slice(work_node,sliced_seq,ExitDipedge,0.5);  // write end dipedge for dipoles

} // SeqElList::slice_node_translate()

void SeqElList::slice_node_default() // slice/translate and add slices to sliced sequence  sliced_seq,  here make always end markers
{
  const element* thick_elem=work_node->p_elem; // work directly on this element  --  to do translation only once, element maybe used several times in sequence


  command_parameter* length_param = return_param_recurse("l",thick_elem); // get original length, value or expression
  if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " length_param=" << length_param << '\n';

  expression* at_expr = work_node->at_expr;
  double at = work_node->at_value;
  if(at_expr==nullptr)  at_expr=expr_from_value_2(at); // make a new expression from the value

  double length = work_node->length;

  expression* l_expr = nullptr;
  if (length_param) l_expr  = length_param->expr;

  int middle=-1;
  if (nslices>1) middle = nslices/2+1; // used to determine after which slice to place the central marker in case of thin slicing

  std::string local_slice_style="collim"; // for passive elements use slice_style "collim"   with slice at start and end

  for (int i=1; i<=nslices; ++i) // loop to place the nslices in the sequence
  {
    element* wirecoll=nullptr;
    if(std::string(thick_elem->base_type->name)=="collimator") wirecoll=create_wire_element(thick_elem,i); // special case wire collimator

    element* slice_i = create_sliced_element(thick_elem,i); // same type, l set to 0 and kept in lrad   (if existing)
    if(wirecoll) for(unsigned int i=0;i<MaTh::WireCollimatorParmList.size();++i) ParameterRemove(MaTh::WireCollimatorParmList[i], slice_i);  // remove wire parameter provided by the wirecoll

    expression* thin_at_expr=nullptr;
    if (fabs(at_shift(nslices,i,local_slice_style))>0.0)
    {
      if( MaTh::iMoreExpressions<1 ) thin_at_expr = compound_expr(at_expr,work_node->at_value,"+",  nullptr, length*at_shift(nslices,i,local_slice_style),1);  // use length and shift values, no expressions
      else thin_at_expr = compound_expr(at_expr,work_node->at_value,"+",scale_expr(l_expr,at_shift(nslices,i,local_slice_style)),length*at_shift(nslices,i,local_slice_style),1); // use length expressions
    }
    else
    {
      if (at_expr) thin_at_expr = clone_expression(at_expr);
    }


    if (i==middle) // create and place new marker with name of thick_elem  in the middle = poition of thick_elem
    {
      element* middle_marker=new_marker_element(thick_elem->name,thick_elem);
      place_node_at(work_node, sliced_seq, middle_marker,at_expr);
    }
    if(slice_i) place_node_at(work_node, sliced_seq, slice_i,thin_at_expr);
    if(wirecoll) place_node_at(work_node, sliced_seq, wirecoll,thin_at_expr);
  }
} // SeqElList::slice_node_default()

void SeqElList::slice_attributes_to_slice(command* cmd,const element* thick_elem)
{ // looks for attributes like kick that need slicing, modify them in cmd used to construct the sliced element
  const std::vector<std::string> attributes_to_slice = { "kick", "hkick", "vkick", "chkick", "cvkick", "current"};
  ElmAttr theElmAttr(thick_elem);
  std::vector<std::string> active_par_to_be_used=theElmAttr.get_list_of_active_attributes();
  for(unsigned i=0;i<active_par_to_be_used.size();++i)
  {
    bool found=false;
    std::string attribute_name=active_par_to_be_used[i];
    for(unsigned int j=0;j<attributes_to_slice.size();++j) // all active attributes
    {
      if(attribute_name==attributes_to_slice[j]) // found in attributes_to_slice
      {
        found=true;
        break;
      }
    }
    if(found)
    {
      name_list* nl=cmd->par_names;
      const int ei=name_list_pos(attribute_name.c_str(),nl);
      if(ei>-1)
      {
        command_parameter* cp=cmd->par->parameters[ei];
        if (cp->type==k_double_array) // double array,  expr_list as used for ilnorm; scale every value, or expression
        {
          if (cp->double_array)
          {
            if (cp->expr_list) // calculate the values
            {
              for (int ei = 0; ei < cp->double_array->curr; ++ei)
              {
                if (ei < cp->expr_list->curr && cp->expr_list->list[ei] != nullptr)
                {
                  cp->expr_list->list[ei] = compound_expr(cp->expr_list->list[ei],0.,"/",nullptr,nslices,1); // scale expression in list
                }
              }
              for (int ei = 0; ei < cp->double_array->curr; ++ei) cp->double_array->a[ei]/=nslices; // scale value in array
            }
          }
        }
        else // single expression/value
        {
          if(cp->expr) cp->expr = compound_expr(cp->expr,0.,"/",nullptr,nslices,1);
          else cp->double_value /= nslices;
        }
      }
    }
  }
}

element* SeqElList::create_sliced_element(const element* thick_elem, int slice_no)
{
  element *sliced_elem_parent;
  if (thick_elem == thick_elem->parent) return nullptr; // no further parent to consider
  else sliced_elem_parent = create_sliced_element(thick_elem->parent,slice_no); // recursively slice parent       --- element specific
  std::string parent_or_base= thick_elem->base_type->name;
  if(sliced_elem_parent) parent_or_base=sliced_elem_parent->name;
  if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << std::left << std::setw(MaTh::el_name_maxlen) << thick_elem->name << " parent_or_base= " << std::setw(MaTh::el_type_maxlen) << parent_or_base << " slice_no=" << slice_no << " nslices=" << nslices << '\n';
  element *sliced_elem;
  if( (sliced_elem = theSliceList->find_slice(thick_elem,slice_no)) )
  {
    return sliced_elem; // already done
  }
  const command_parameter* at_param = return_param_recurse("at",thick_elem);  // handle parent with possibly different slice number than child
  if(!at_param && thick_elem == thick_elem->parent) slice_no=nslices=1; // do not slice this one
  if(slice_no > nslices && thick_elem!=thick_elem->parent ) slice_no=1; // check, but not for base classes
  if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << std::left << std::setw(MaTh::el_name_maxlen) << thick_elem->name << " " << std::setw(MaTh::el_type_maxlen) << work_node->base_name << " slice_no=" << slice_no << " nslices=" << nslices << " is parent=" << (thick_elem==thick_elem->parent) << '\n';
  const command_parameter* length_param  = return_param_recurse("l",thick_elem);

  command* cmd = new_cmdptr( thick_elem );
  copy_params_from_elem(cmd,thick_elem,str_v_join(MaTh::DoNotCopy,{"kill_ent_fringe","kill_exi_fringe"}));
  set_lrad(cmd,length_param,nslices); // keep l  as lrad
  if(nslices>1) slice_attributes_to_slice(cmd,thick_elem); // like kick
  finish_make_sliced_elem(sliced_elem, thick_elem, cmd, parent_or_base, slice_no);
  ParameterRemove("l", sliced_elem);
  return sliced_elem;
}

void SeqElList::finish_make_sliced_elem(element*& sliced_elem, const element* thick_elem,command* cmd,const std::string parent_name,int slice_no)
{ // final steps in making sliced element, set sliced name, make element from cmd and put to theSliceList
  std::string thick_elem_name=thick_elem->name;
  if(thick_elem_name[0]==MaTh::ExtraChar)
  {
    thick_elem_name=thick_elem_name.substr(1);
  }
  std::string thin_name;
  if ( MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " thick_elem->length=" << thick_elem->length << std::endl;
  if (nslices==1 && slice_no==1) thin_name=thick_elem_name;
  else thin_name = make_thin_name(thick_elem_name,slice_no); // add slice number
  sliced_elem = my_El_List->my_make_element(thin_name,parent_name,cmd,-1);
  theSliceList->put_slice(thick_elem,sliced_elem);
}

node* SeqElList::copy_thin(node* work_node) // this copies an element node and sets the length to zero and lrad to the length to be used for "copying" optically neutral elements
{
  if ( MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << "  " << std::setw(MaTh::par_name_maxlen) << work_node->name << " " << std::setw(MaTh::el_type_maxlen) << work_node->base_name << " thin_node->length=" << work_node->length << " l=" << el_par_value("l",work_node->p_elem) << std::endl;
  node* thin_node = nullptr;
  thin_node = clone_node(work_node, 0);
  if (el_par_value("l",work_node->p_elem)>zero)
  {
    if ( MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << "  " << std::setw(MaTh::par_name_maxlen) << work_node->name << " had length, remove" << '\n';
    thin_node->p_elem = create_sliced_element(work_node->p_elem,1); // keep original length as lrad
  }
  thin_node->length=0;
  thin_node->p_elem->length=0;
  return thin_node;
}

void SeqElList::slice_node() // main steering, decides how to split an individual node and sends it onto the sliced_seq builder
{
  static const std::vector<std::string> element_to_slice_and_tanslate = { "elseparator", "octupole", "quadrupole", "rbend", "sbend", "sextupole", "solenoid" };
  const element* thick_elem=work_node->p_elem;
  nslices = get_slices_from_elem(thick_elem);

  if(verbose>1) std::cout << '\n'  << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " slice_style=\"" << slice_style << "\"" << " nslices=" << nslices;
  if (thick_elem) // look at the element of this node to see what to do with slicing
  {
    if(verbose>1) std::cout << " now see what to do with work_node=" << std::left << std::setw(MaTh::par_name_maxlen) << work_node->name <<  " depending on its base=" << std::setw(MaTh::par_name_maxlen) << work_node->base_name << std::right << '\n';
    const double eps=1.e-15; // used to check if a value is compatible with zero
    bool IsWireCollimator = ( strcmp(work_node->base_name,"collimator") == 0 && return_param_recurse("current",thick_elem) );
    if ( fabs(el_par_value("l",thick_elem)) <eps ) // if the length is compatible with zero copy it directly
    {
      if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " zero length place directly copy of element " << work_node->name << '\n';
      add_node_at_end_of_sequence( copy_thin(work_node),sliced_seq );
    }
    else if( nslices==0)
    {
      if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " number of slices is 0, place as is " << work_node->name << '\n';
      add_node_at_end_of_sequence( work_node,sliced_seq );
    }
    else if ( strcmp(work_node->base_name,"matrix") == 0 ) add_node_at_end_of_sequence(work_node,sliced_seq); // Take matrix as it is, including any length
    else if ( strcmp(work_node->base_name,"drift")  == 0 ) ; // do nothing for drift
    // now the element types that can be sliced
    else if (NameIsInList(work_node->base_name,element_to_slice_and_tanslate))
    {
      slice_node_translate(); // slice/translate and add slices to sliced sequence
    }
    else if(nslices>1 || IsWireCollimator) // in case of wire collimator translate also single slice
    {
      slice_node_default(); // and add slices to sliced sequence
    }
    else // single slice
    {
      node* thin_node;
      if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " nslices=" << nslices << " place copy of work_node " << work_node->name << '\n';
      if(MaTh::iMakeEndMarkers) place_start_or_end_marker(true); // start
      add_node_at_end_of_sequence(thin_node = copy_thin(work_node),sliced_seq);
      if (strcmp(thin_node->p_elem->base_type->name, "rfcavity") == 0 &&
          find_element(thin_node->p_elem->name, sliced_seq->cavities) == nullptr)
        add_to_el_list(&thin_node->p_elem, 0, sliced_seq->cavities, 0);     // special cavity
      if (strcmp(thin_node->p_elem->base_type->name, "crabcavity") == 0 &&
          find_element(thin_node->p_elem->name, sliced_seq->crabcavities) == nullptr)
        add_to_el_list(&thin_node->p_elem, 0, sliced_seq->crabcavities, 0); // special crab cavity
      if(MaTh::iMakeEndMarkers) place_start_or_end_marker(false); // end
    }
  } // done with case where thick_elem  is defined
  else if (work_node->p_sequ) // nested sequence (not flattened), slice and add the subsequence, this case is checked in fivecell.madx
  {
    sequence* sub_thin=theSequenceList->slice_sequence(slice_style,work_node->p_sequ);
    node* sub_node = new_sequ_node(sub_thin, work_node->occ_cnt);
    sub_node->length = 0;
    sub_node->at_value = work_node->at_value;
    if (sub_node->at_expr) sub_node->at_expr = clone_expression(work_node->at_expr);
    if(verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " place the sliced sub-sequence " << work_node->p_sequ->name << '\n';
    add_node_at_end_of_sequence(sub_node,sliced_seq);
  }
  else fatal_error("node is not element or sequence",work_node->base_name); // completely unknown, error
}

void SeqElList::place_thick_slice(const element* thick_elem, element* sliced_elem, const int i) // make nodes for the _s, _b  pieces  and place them in the sequence
{
  if(sliced_elem==nullptr) return; // nothing to place
  const int n=nslices-1; // in case of thick slices,
  SliceDistPos SP(n, slice_style==std::string("teapot") );
  if(verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " " << thick_elem->base_type->name << " " << sliced_elem->name << " nslices=" << nslices << " n=" << n << " start from thick_elem " << thick_elem->name << " MaTh::iMoreExpressions=" << MaTh::iMoreExpressions << std::endl;

  expression* l_par_expr=my_get_param_expression(thick_elem, "l"); // with this l_par_expr should not be NULL
  expression* at_expr = clone_expression(work_node->at_expr);
  double at = work_node->at_value;

  if(verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " at_expr " << my_dump_expression(at_expr) << '\n';

  double rel_shift;

  if(MaTh::iMoreExpressions<2)
  { // use double value for shift, no expression
    if(nslices==1)       rel_shift=0; // single thick piece remains in centre
    else if(i==1)        rel_shift=-0.5 + SP.delta/2.; // entry
    else if(i==nslices)  rel_shift= 0.5 - SP.delta/2.; // exit
    else                 rel_shift=-0.5 + SP.delta + (i-1.5)*SP.Delta; // body
    if(verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " rel_shift=" << rel_shift << '\n';
    if(MaTh::iMoreExpressions<1) at_expr = compound_expr(at_expr,at, "+", nullptr, el_par_value("l",thick_elem) *rel_shift,1 );  // use length and shift values, no expressions
    else at_expr = compound_expr(at_expr,at, "+", scale_expr(l_par_expr,rel_shift),  0 ,1); // MaTh::iMoreExpressions==1, use length expression and shift value
  }
  else
  { //
    expression* rel_shift_expr;
    if(nslices==1)
    {
      rel_shift_expr = new_expression("0", nullptr); // single thick piece remains in centre
    }
    else if(i==1) // entry piece -1/2 + 1./(2.*SP.delta_inv)
    {
      std::string tstr1="-1/2";
      rel_shift_expr = compound_expr(new_expression(tstr1.c_str(),nullptr), 0., "+", new_expression(SP.delta_half_str.c_str(),nullptr), 0,1); // entry
    }
    else if(i==nslices) // 0.5 - 1./(2.*SP.delta_inv); // exit
    {
      std::string tstr1="1/2";
      rel_shift_expr = compound_expr(new_expression(tstr1.c_str(),nullptr), 0., "-", new_expression(SP.delta_half_str.c_str(),nullptr), 0,1); // exit
    }
    else // -0.5 + 1./SP.delta_inv + (i-1.5)*SP.Delta; // body
    {
      std::string tstr1="-1/2";
      std::string tstr2=SP.delta_str+"+"+std::to_string(2*i-3)+"*"+SP.Delta_half_str;
      rel_shift_expr = compound_expr(new_expression(tstr1.c_str(),nullptr), 0., "+", new_expression(tstr2.c_str(),nullptr), 0,1); // entry
    }
    // multiply with lpar
    rel_shift_expr=compound_expr(l_par_expr,at, "*", rel_shift_expr,  0,1 );
    if(verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " rel_shift_expr " << my_dump_expression(rel_shift_expr) << '\n';
    at_expr = compound_expr(at_expr,at, "+", rel_shift_expr, 0, 1); // this also updates the value
  }

  place_node_at(work_node,sliced_seq,sliced_elem,at_expr);
}

void  SeqElList::place_start_or_end_marker(const bool at_start)
{
  const element* thick_elem=work_node->p_elem;
  if(verbose>1) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " work_node " << work_node->name << " at_start=" << at_start << '\n';
  std::string AddToName;
  double rel_shift;
  if(at_start)
  {
    AddToName="_mken"; // for marker at entry
    rel_shift=-0.5;    // -0.5 length from centre
  }
  else
  {
    AddToName="_mkex"; // for marker at exit
    rel_shift= 0.5;    // +0.5 length from centre
  }
  std::string thick_elem_name=thick_elem->name;
  if(thick_elem_name[0]==MaTh::ExtraChar)
  {
    thick_elem_name=thick_elem_name.substr(1); // without the ExtraChar
  }
  element* start_end_marker=new_marker_element(thick_elem_name+AddToName,thick_elem);
  place_thin_slice(work_node,sliced_seq,start_end_marker,rel_shift);
}

//--------  SequenceList

SequenceList::SequenceList()  // constructor
{
  if(theSliceList   ==nullptr) theSliceList    = new ElementListWithSlices(MaTh::Verbose);
  if(theRbendList   ==nullptr) theRbendList    = new ElementListWithSlices(MaTh::Verbose);
  if(theBendEdgeList==nullptr) theBendEdgeList = new ElementListWithSlices(MaTh::Verbose);
  if(my_El_List==nullptr)      my_El_List      = new my_Element_List();
}

sequence* SequenceList::find_sequ(sequence* thick_sequ) // check if thick_sequ is already in my_sequ_list_vec
{
  if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " my_sequ_list_vec.size()=" << my_sequ_list_vec.size() << '\n';
  for(unsigned int i=0; i<my_sequ_list_vec.size(); ++i)
  {
    if ( my_sequ_list_vec[i] == thick_sequ )
    {
      return thick_sequ;
    }
  }
  return nullptr;
}

void SequenceList::put_sequ(sequence* thick_sequ)
{
  my_sequ_list_vec.push_back(thick_sequ);
  if(MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " my_sequ_list_vec.size()=" << my_sequ_list_vec.size() << '\n';
  return;
}

void SequenceList::Print(std::ostream& StrOut) const
{
  StrOut << "SequenceList::Print() currently " << my_sequ_list_vec.size() << " defined:" << '\n';
  for(unsigned int i=0; i<my_sequ_list_vec.size(); ++i) StrOut << " " << my_sequ_list_vec[i]->name;
  StrOut << '\n';
  return;
}

void SequenceList::Reset()
{
  if ( MaTh::Verbose>1) std::cout << __FILE__ << " " << __PRETTY_FUNCTION__ << " line " << std::setw(4) << __LINE__ << " before reset my_sequ_list_vec.size()=" << my_sequ_list_vec.size() << '\n'; // debug
  my_sequ_list_vec.resize(0);

  delete theSliceList;
  delete theRbendList;
  delete theBendEdgeList;
  delete my_El_List;

  theSliceList    = new ElementListWithSlices(MaTh::Verbose);
  theRbendList    = new ElementListWithSlices(MaTh::Verbose);
  theBendEdgeList = new ElementListWithSlices(MaTh::Verbose);
  my_El_List      = new my_Element_List();

}

sequence* SequenceList::slice_sequence(const std::string slice_style,sequence* thick_sequ,const std::string LastSequenceSliced,const std::string LastStyle) // make recursively a sliced sequence out of the thick_seque
{
  static std::string thick_sequ_name_last;;
  std::string thick_sequ_name("null");
  if(thick_sequ) thick_sequ_name=thick_sequ->name;

  if(thick_sequ_name==LastSequenceSliced && slice_style!=LastStyle)
  {
    Reset();
  }

  sequence* sliced_seq=find_sequ(thick_sequ); // check if already in list

  if(sliced_seq)
  {
    if(MaTh::Verbose) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " sequence " << thick_sequ->name << " was already sliced" << std::endl;
    if(MaTh::Verbose>1) Print();
    return sliced_seq; // do nothing if the sequence was already sliced
  }

  const std::string name = thick_sequ->name;
  std::cout << "makethin: slicing sequence : " << name << '\n';
  if(MaTh::Verbose>1)
  {
    int level=1;
    if(MaTh::Verbose>2) level=2;
    std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " my_dump_sequence thick_sequ " << my_dump_sequence(thick_sequ,level) << std::endl; // dump level 2, without nodes/elements
  }

  sliced_seq = new_sequence(name.c_str(), thick_sequ->ref_flag);
  sliced_seq->start = nullptr;
  sliced_seq->share = thick_sequ->share;
  sliced_seq->nested = thick_sequ->nested;
  sliced_seq->length = sequence_length(thick_sequ);
  sliced_seq->refpos = permbuff(thick_sequ->refpos);
  sliced_seq->beam = thick_sequ->beam;
  if (sliced_seq->cavities != nullptr)  sliced_seq->cavities->curr = 0;
  else sliced_seq->cavities = new_el_list(100);
  if (sliced_seq->crabcavities != nullptr)  sliced_seq->crabcavities->curr = 0;
  else sliced_seq->crabcavities = new_el_list(100);

  SeqElList theSeqElList(name, slice_style, sliced_seq, thick_sequ->start,this);
  while(theSeqElList.current_node() != nullptr) // in current sequence, loop over nodes
  {
    theSeqElList.slice_node(); // decides what to do with current node :  slice_node_translate, slice_node_default or nothing
    if (theSeqElList.current_node() == thick_sequ->end)
    {
      break;
    }
    if(theSeqElList.current_node()->p_elem!=nullptr){
      if(strcmp(theSeqElList.current_node()->p_elem->base_type->name, "rfcavity")==0 &&
        find_element(theSeqElList.current_node()->p_elem->name, sliced_seq->cavities) == nullptr){
        add_to_el_list(&theSeqElList.current_node()->p_elem, 0, sliced_seq->cavities, 0);
      }
    }
    theSeqElList.current_node(theSeqElList.current_node()->next); // set current_node
  }
  sliced_seq->end->next = sliced_seq->start;


  int pos=0;
  if ((pos = name_list_pos(name.c_str(), sequences->list)) < 0) // move the pointer in the sequences list to point to our thin sequence
  {
    fatal_error("unknown sequence sliced:", name.c_str());
  }
  else
  {
    sequences->sequs[pos]= sliced_seq; // pointer moved ok, delete_sequence(thick_sequ)
  }

  thick_sequ_name_last=thick_sequ_name;

  put_sequ(thick_sequ); // Slicing done for this sequence. Add to list of sequences sliced
  if(MaTh::Verbose) std::cout << __FILE__ << " " << __FUNCTION__ << " line " << std::setw(4) << __LINE__ << " before print theSeqElList" << std::endl;
  if(MaTh::Verbose) theSeqElList.Print(); // print final list

  return sliced_seq;
} // slice_sequence
