/*
 mad_mkthin.cpp

 Thick to thin lens converter makethin. Helmut Burkhardt

 Major steps
 2005 : Standard selection SELECT,FLAG=makethin,RANGE=range,CLASS=class,PATTERN=pattern[,FULL][,CLEAR];
 Implementation of slicing for solenoids
 2012 : Extension of TEAPOT slicing to n>4
 2013 : Keep thick elements if slice number <1, Code now C++,
        Thick slicing for quadrupoles
	    Automatic generation of dipedge elements for dipoles
 2014 : Thick bend slicing, with or without dipedge

 here test version with
 static const bool dipedge_h1_h2_fl=true;  // normally false to avoid potentially non-simplectic partial higher order in dipedge. Optionally true as requested by Andrea Latina in 10/2014
 static const bool kill_fringe_fl=true;    // requested by Laurent et al., probably redundant, should be sufficient to check existance of non-default h1,e1; h2,e2 parameters
 static const bool rbend_to_sbend_fl=true; // translate all rbends to sbends, even if  makedipedge=false

 Early versions in 2001, 2002 by Mark Hayes

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

static const char EOL='\n'; // to avoid problems with older C++ compilers, where character literals are integer

//------------------------------- forward declarations --------------

class SliceDistPos // defines the distances and positions, depending on number of slices and slicing style
{
public:
  SliceDistPos(const int n,const bool teapot_fl); // constructor
  ~SliceDistPos() {}; // (empty) destructor
  void Print() const;
  double delta;
  double Delta;
private:
  int n; // number of slices
  bool teapot_fl;
};

class OneElementWithSlices // One Element with Slices       used to work on slices, derived from thick_elem which is not modified - declared as constant
{
public:
  OneElementWithSlices(const element* thick_elem,element* thin_elem); // constructor
  ~OneElementWithSlices() {}; // (empty) destructor
  const element* thick_elem; // pointer to the thick element
  std::vector<element*> sliced_elem; // pointer(s) to the one or several slices
};

class ElementListWithSlices
{
public:
  std::vector<OneElementWithSlices*> VecElemWithSlices; // vector of thick elements+slices
  ElementListWithSlices(unsigned int verbose); // constructor
  ~ElementListWithSlices(); // destructor
  void put_slice(const element* thick_elem, element* thin_elem); // add thin_elem to VecElemWithSlices
  element* find_slice(const element* thick_elem, const int slice); // find address of thin slice by slice number for thick_elem
  element* find_slice(const element* thick_elem, const std::string& name); // find address of thin slice by slice name for thick_elem
  void Print() const;
  void PrintCounter() const;
private:
  int find_thick(const element* thick_elem); // find thick_element in VecElemWithSlices, <0 means not found,  used inside find_slice
  unsigned int verbose;
  unsigned int get_thin_calls,get_thin_iteractions; // to monitor the search (in-) efficiency find_slice
  int ilast1,ilast2; // keep last two found find_slice, useful in recursive searches which switch between slices and parents
};

class SeqElList // sequence with elements considered for slicing
{
public:
  SeqElList(const std::string& seqname,const std::string& thin_style,sequence* thick_sequ,sequence* thin_sequ,node* thick_node); // constructor
  ~SeqElList(); // destructor
  void Print() const;
  void slice_node();
  node* thick_node; // current node
private:
  double  simple_at_shift(const int slices, const int slice_no);
  double  teapot_at_shift(const int slices, const int slice_no);
  double  collim_at_shift(const int slices, const int slice_no);
  double default_at_shift(const int slices, const int slice_no);
  double at_shift(const int slices, const int slice_no,const std::string& local_thin_style); // return at relative shifts from centre of unsliced magnet
  int translate_k(command_parameter* *kparam, command_parameter* *ksparam,const command_parameter* angle_param, command_parameter* kn_param, command_parameter* ks_param);
  element* create_thick_slice(element* thick_elem,const int slice_no);
  element* create_sliced_magnet(const element* thick_elem, int slice_no,bool ThickSLice);
  element* create_thin_solenoid(const element* thick_elem, int slice_no);
  element* create_thin_elseparator(const element* thick_elem, int slice_no);
  element* create_thin_obj(const element* thick_elem, int slice_no);
  void slice_this_node(); // called in loop over nodes which can be sliced, makes slices, adds them to the sliced sequence
  node* copy_thin(node* thick_node);
  ElementListWithSlices *theSliceList, *theBendEdgeList; // Elements, of various types; consider to make separate lists for quadrupoles, sext ..  to speed up search
  std::string seqname; // name of the sequence
  std::string thin_style;
  sequence *thick_sequ, *thin_sequ;
  unsigned int verbose;
  const double eps;
  bool MakeDipedge; // translate dipoles   to    dipedge, dipole without edge effects, dipedge
};

class SequenceList
{
public:
  void put_sequ(sequence* thick_sequ); // add a sequence to the sequence list
  sequence* get_sequ(sequence* thick_sequ); // check if thick_sequ is already there, if yes return the pointer to it,  used to check if the sequence was already sliced
private:
  std::vector<sequence*> my_sequ_list_vec; // list of sequences
};

//------------------------------- source code --------------

static const int        k_logical=0;
static const int            k_int=1;
static const int         k_double=2;
static const int        k_cstring=3;
// static const int        k_min_max=4;
static const int     k_int_array=11;
static const int  k_double_array=12;
static const int k_cstring_array=13;
using namespace std;

static const bool dipedge_h1_h2_fl=true;  // normally false to avoid potentially non-simplectic partial higher order in dipedge. Optionally true as requested by Andrea Latina in 10/2014
static const bool kill_fringe_fl=true;    // requested by Laurent et al., probably redundant, should be sufficient to check existance of non-default h1,e1; h2,e2 parameters
static const bool rbend_to_sbend_fl=true; // translate all rbends to sbends, even if  makedipedge=false

// check general options
inline bool debug_fl()    { return get_option( const_cast<char*>("debug") ); }
inline bool verbose_fl()  { return get_option( const_cast<char*>("verbose") ); }
inline bool thin_foc_fl() { return get_option( const_cast<char*>("thin_foc") ); }
inline bool rbarc_fl()    { return get_option( const_cast<char*>("rbarc") ); } // by default on, then use (reduced) length of rbends

static char* CopToNewC_String(const string& s) // copy chars of string to a new c-string Stroustrup3 p.590
{
  // char* p=new char[s.length()+1]; // note: +1  (for \0)
  char* p = static_cast<char*>(mymalloc_atomic("CopToNewC_String",s.length()+1)); // allocation for C, allowing for myfree
  memcpy(p,s.c_str(),s.length()+1);   // does not copy the terminator \0
  return p;
}


//--------------------  const versions of mad utility routines ---- start

/* returns parameter if it has been modified, otherwise NULL */
static const command_parameter*
return_param(const char* par,const element* elem) //hbu const version
{
  int index;
  /* don't return base type definitions */
  if (elem==elem->parent) return NULL;

  if ((index = name_list_pos(par,elem->def->par_names))>-1
      && elem->def->par_names->inform[index] > 0)
    return elem->def->par->parameters[index];
  return NULL;
}

/* returns parameter if it has been modified, otherwise NULL  - recursively */
static const command_parameter*
return_param_recurse(const char* par,const element* elem) //hbu const version
{
  const command_parameter* param = return_param(par,elem);
  if (param) return param;
  if (elem!=elem->parent)
    return return_param_recurse(par,elem->parent);
  return NULL;
}

static command_parameter*
clone_command_parameter(const command_parameter* p) //hbu const version, also protect that does not crash if NULL
{
  if(p==NULL) return NULL; //hbu protect from crashing
  command_parameter* clone = new_command_parameter(const_cast<char*>(p->name), p->type); //hbu
  clone->call_def = p->call_def;
  switch (p->type)
  {
    case 4:
      clone->c_min = p->c_min;
      clone->c_max = p->c_max;
      clone->min_expr = clone_expression(p->min_expr);
      clone->max_expr = clone_expression(p->max_expr);
    case 0:
    case 1:
    case 2:
      clone->double_value = p->double_value;
      clone->expr = clone_expression(p->expr);
      break;
    case 3:
      clone->string = p->string;
      clone->expr = NULL;
      break;
    case 11:
    case 12:
      clone->double_array = clone_double_array(p->double_array);
      clone->expr_list = clone_expr_list(p->expr_list);
      break;
    case 13:
      clone->m_string = clone_char_p_array(p->m_string);
  }
  return clone;
}

static void
add_cmd_parameter_clone(command* cmd,const command_parameter *param,const char* par_name,const int inf) // hbu add an identical copy (clone) of param to cmd,  here const version
{
  if(param)
  {
    cmd->par->parameters[cmd->par->curr] = clone_command_parameter(param); /* set current to identical copy (clone) of param */
    add_to_name_list(const_cast<char*>(par_name),inf,cmd->par_names); //hbu const_cast<char*>(
    cmd->par->curr++;
  }
}

static double
command_par_special(const char* parameter,const element* el) //hbu const version
/* construct missing tilt from normal and skew  */
{
  double val = zero;

  if (strcmp(parameter, "tilt") == 0)
  {
    if ((val = command_par_value("tilt", el->def)) == zero)
    {
      val = zero;
    }
  }
  else val = command_par_value(parameter, el->def);
  return val;
}

static int
element_vector(const element* el,const char* par,double* vector) //hbu const version
/* returns length + vector of parameter par for element el */
{
  int i, l = 0;
  double_array* da;
  expr_list* ell;
  if ((i = name_list_pos(par, el->def->par_names)) > -1)
  {
    if ((da = el->def->par->parameters[i]->double_array) != NULL)
    {
      if ((ell = el->def->par->parameters[i]->expr_list) != NULL)
        update_vector(ell, da);
      l = da->curr;
      copy_double(da->a, vector, l);
    }
  }
  return l;
}

static double
el_par_value(const char* par,const element* el) //hbu const version
/* returns an element parameter value */
{
  int k = 0; // , n; not used
  char tmp[8];
  double val = zero, angle = zero, l, vec[100];
  double fact = strcmp(el->base_type->name, "rbend") == 0 ? one : zero;
  int mult = strcmp(el->base_type->name, "multipole") == 0 ? 1 : 0;
  int mark = strcmp(el->base_type->name, "marker") == 0 ? 1 : 0;
  if (fact != zero || strcmp(el->base_type->name, "sbend") == 0) /* bend */
  {
    if ((l = command_par_value("l", el->def)) == zero)
      fatal_error("bend with zero length:",el->name);
    angle = command_par_value("angle", el->def);
    if (strcmp(par, "angle") == 0)  val = angle;
    else if (strcmp(par, "tilt") == 0)
      val = command_par_value("tilt", el->def);
    else if (strcmp(par, "k0") == 0) val = command_par_value("k0", el->def);
    else if (strcmp(par, "k0s") == 0) val = command_par_value("k0s", el->def);
    else if (strcmp(par, "l") == 0)
    {
      if (fact != zero && get_option( const_cast<char*>("rbarc") ) && angle != zero) //hbu const_cast
        val = l * angle / (two * sin(angle/two));
      else val = l;
    }
    else if (strcmp(par, "e1") == 0)
      val = command_par_value("e1", el->def);/* + fact * angle / two; dipole_bv kill initiative SF TR FS */
    else if (strcmp(par, "e2") == 0)
      val = command_par_value("e2", el->def);/* + fact * angle / two; dipole_bv kill initiative SF TR FS */
    else if (strcmp(par, "rhoinv") == 0) val = angle / l;
    else if (strcmp(par, "blen") == 0) val = l;
    else val = command_par_value(par, el->def);
  }
  /* all elements except bends */
  else if (strcmp(par, "rhoinv") == 0) val = zero;
  else if (strcmp(par, "blen") == 0) val = zero;
  else if (mark) /* marker */
  {
    if ((l = command_par_value("l", el->def)) != zero)
      fatal_error("marker with nonzero length:",el->name);
    val = command_par_value(par, el->def);
  }
  else if (mult)  /* multipole */
  {
    if (strcmp(par, "l") == 0) val = zero;
    else if (par[0] == 'k' && isdigit(par[1]) && par[strlen(par)-1] == 'l')
    /* single component requested for multipole */
    {
      if (strchr(par, 's')) strcpy(tmp, "ksl");
      else                  strcpy(tmp, "knl");
      sscanf(&par[1], "%d", &k);
      if (element_vector(el, tmp, vec) > k)  val = vec[k]; // n = not used
    }
    else val = command_par_value(par, el->def);
  }
  else val = command_par_special(par, el);
  /* extra code for kickers */
  if (val == zero && strcmp(el->base_type->name, "hkicker") == 0)
  {
    if (strcmp(par,"hkick") == 0)
      val = command_par_value("kick", el->def);
    else if (strcmp(par,"kick") == 0)
      val = command_par_value("hkick", el->def);
  }
  else if (val == zero && strcmp(el->base_type->name, "vkicker") == 0)
  {
    if (strcmp(par,"vkick") == 0)
      val = command_par_value("kick", el->def);
    else if (strcmp(par,"kick") == 0)
      val = command_par_value("vkick", el->def);
  }
  return val;
}

//hbu add_cmd_parameter_new written by me and only used here, no need to mad_cmdpar.c,h version not needed
static void
add_cmd_parameter_new(command* cmd,const double par_value,const char* par_name,const int inf) //hbu const version
{
  cmd->par->parameters[cmd->par->curr] = new_command_parameter(const_cast<char*>(par_name), 2);
  cmd->par->parameters[cmd->par->curr]->double_value = par_value;
  add_to_name_list(const_cast<char*>(par_name),inf,cmd->par_names);
  cmd->par->curr++;
}

//--------------------  const versions of mad utility routines ---- end

static bool NameIsInList(const string& Name,const vector<string>& List)
{
  bool found=false;
  for(unsigned int i=0;i < List.size(); ++i)
  {
    if(Name==List[i])
    {
      found=true;
      break;
    }
  }
  return found;
}

static string my_dump_expression(expression* ex) // dump_expression in mad_expr.c only prints the value, here show also the expression and check for NULL
{
  ostringstream ostr;
  ostr << setprecision(15) << "expression ";
  if(ex==NULL) ostr << " is NULL";
  else
  {
    if(ex->string) ostr << " string=" << left << setw(20) << ex->string << right;
    ex->value = expression_value(ex, 2);
    ostr << " value=" << ex->value;
  }
  return ostr.str();
}

static double my_get_int_or_double_value(const element* el,const char* parnam,bool &found) // works for integer and double, also useful as   my_get_int_or_double_value(el->base_type,char* parnam);  to get the default values
{
  // just returning  el_par_value(parnam,base_el);   or  el_par_value_recurse(parnam,base_el);
  // is not good enough, gets 0 for integer parameters and in some cases too specific - like checks for non-zero length of dipoles  which does not allow to compare with the base type for dipoles
  // rather do all directly here, descending from el / el->def / el->def->par to parameters[i], loop through them and look at integer and double values
  // in case of expression uses the expression_value
  double val=0;
  found=false;
  unsigned int verbose=0;
  if(verbose_fl()) verbose=2;
  // verbose=3;
  if(el && el->def && el->def->par)
  {
    command_parameter_list* pl=el->def->par;
    if(verbose>2) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " el->name=" << setw(15) << left << el->name;
    for (int i = 0; i < pl->curr; ++i)
    {
      if(pl->parameters[i])
      {
        command_parameter* cp=pl->parameters[i];
        if( string(cp->name) == string(parnam) )
        {
          if(cp->expr)
          {
            val = expression_value(cp->expr,2);
            if(verbose>2) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " el->name=" << setw(15) << left << el->name << " use the expression_value=" << val << EOL;
            found=true;
          }
          else switch (cp->type)
          {
            case k_int:    //    int value of expression, actually same as double value
              if(verbose>2) cout << "    int ";
              found=true;
              val = cp->double_value;
              break;
            case k_double: // double value of expression
              if(verbose>2) cout << " double ";
              found=true;
              val = cp->double_value;
              break;
          }
        }
      }
    }
  }
  if(verbose>2)
  {
    cout << " parameter " << setw(15) << parnam;
    if(found) cout << "     found"; else      cout << "         not found";
    cout << right << " val=" << val << EOL;
  }
  return val;
}

static string my_dump_name_list(const name_list* nl) // name_list defined in mad_name.h
{
  ostringstream ostr;
  ostr << setprecision(15) << "my_dump_name_list  name          i   j   inform    max=" << setw(2) << nl->max << " curr=" << setw(2) << nl->curr << " name=" << setw(30) << nl->name << " list address=" << nl << endl; // max typically defined in new_name_list in mad_name.c
  if(nl==NULL) ostr << " is NULL";
  else
  {
	for (int i = 0; i < nl->curr; ++i)
	{
	  const int j=nl->index[i];
	  string nl_name="NULL";
	  if(nl->names[j]) nl_name=string(nl->names[j]); else ostr << " *** warning ***  name for " << i << " is NULL, this should not happen" << EOL;
	  if(nl_name.length()>30) ostr << " *** warning ***  name for " << i << " is long, length=" << nl_name.length() << " looks like something was overwritten,  name within quotes =\"" << nl_name << "\"" << EOL;
	  ostr
	  << setw(30) << left << nl_name << right
	  << setw(4) << i
	  << setw(4) << j
	  << setw(4) << nl->inform[j]
	  << EOL;
	}
  }
  return ostr.str();
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

static string my_dump_command_parameter(const command_parameter* cp) // dump_command_parameter in mad_cmdpar.c only prints the value, here show also the expression by calling my_dump_expression and check for NULL
{
  ostringstream ostr;
  ostr << setprecision(15) << "my_dump_command_parameter ";
  if(cp==NULL) ostr << " is NULL";
  else
  {
	ostr << "parameter:" << left << setw(15) ;
	if(cp->name)  ostr << cp->name; else ostr << "name=NULL ";
	ostr << right << " cp->type=" << setw(2) << cp->type;
	ostr << " stamp=" << cp->stamp << " ";
    double default_val=0;
	const double eps=1.e-15; // used to check if strength is compatible with zero
	switch (cp->type)
	{
	  case k_logical:
		ostr << "logical: ";
		if( (int) cp->double_value) ostr << "true"; else ostr << "false";
		ostr << EOL;
		break;
	  case k_int:    // int    value of expression, actually same as double value
	  case k_double: // double value of expression
		if(cp->expr) ostr << my_dump_expression(cp->expr); else ostr << " expression=NULL ";
        if(cp->call_def) default_val=cp->call_def->double_value;
 		if(cp->expr==NULL && fabs(cp->double_value-default_val)>eps)
          ostr << " value=" << setw(10) << cp->double_value << setw(10) << default_val;
		ostr << EOL;
		break;
	  case k_int_array:    // int array,     expr_list
	  case k_double_array: // double array,  expr_list, used for example for Aperture, http://mad.web.cern.ch/mad/madx.old/Introduction/aperture.html
		if (cp->double_array != NULL)
		{
		  if (cp->expr_list != NULL) // calculate the values
		  {
			ostr << "array of " << cp->double_array->curr << "  ";
			for (int ei = 0; ei < cp->double_array->curr; ++ei)
			{
			  if (ei < cp->expr_list->curr && cp->expr_list->list[ei] != NULL)
			  {
				ostr << right << setw(3) << ei << " :" << left << my_dump_expression(cp->expr_list->list[ei]) << right; // show expression and value
			  }
			}
		  }
		}
		ostr << EOL;
		break;
	  case k_cstring:
		ostr << "cstring:";
		if(cp->string) ostr << cp->string; else ostr << " NULL";
		ostr << EOL;
		break;
	  case k_cstring_array: // string array
		dump_char_p_array(cp->m_string);
	  case '?':
		ostr << " cp->type=" << cp->type << " no info dump implemented so far" << EOL;
	}
  }
  return ostr.str();
}

static string my_dump_command_parameter_list(command_parameter_list* pl)
{
  ostringstream ostr;
  ostr << setprecision(15) << "my_dump_command_parameter_list";
  if(pl==NULL) ostr << " is NULL";
  else
  {
	if(pl->name) ostr << " name=" << pl->name; else ostr << " name=NULL";
	ostr << " curr=" << pl->curr << " max=" << pl->max << EOL;
	if(pl->curr > pl->max)
	{
	  ostr << "*** error *** seen in my_dump_command_parameter_list max=" << pl->curr << " > " << " curr" << pl->curr << " set curr back to max" << EOL;
	  pl->curr = pl->max;
	}
	for (int i = 0; i < pl->curr; ++i)
	{
	  ostr << setw(2) << i << " : ";
	  if(pl->parameters[i]) ostr << my_dump_command_parameter(pl->parameters[i]); else ostr << " NULL ";
	}
  }
  return ostr.str();
}

static void SetParameterValue(const char* parnam,element* el,const double val,const int type=k_double) // set value and type, by default double
{
  const int ei=name_list_pos(parnam,el->def->par_names);
  if(ei > -1)
  {
    command_parameter* cp=el->def->par->parameters[ei];
    if(cp)
    {
      if(verbose_fl())
      {
        cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " parameter " << parnam
        << " was double_value=" << cp->double_value
        << " and type=" << cp->type;
        if(cp->expr) cout << " has " << my_dump_expression(cp->expr); else cout << " no expression";
        cout << " set to val=" << val
             << " and type=" << type << EOL;
      }
      if(cp->expr) cp->expr=NULL; // remove any expression
      cp->double_value=val; // set the double value
      cp->type=type; // set the type value
    }
  }
  else cout << " *** warning *** SetParameterValue for parameter " << parnam << " failed for " << el->name << " parameter not in element name_list " << EOL;
}

static void ParameterTurnOn(const char* parnam,element* el) // request that this parameter is written to output
{
  const int ei=name_list_pos(parnam,el->def->par_names);
  if(ei > -1) el->def->par_names->inform[ei]=1; // Turn on by setting inform to 1
  else cout << " *** warning *** ParameterTurnOn for parameter " << parnam << " failed for " << el->name << " parameter not in element name_list " << EOL;
}

static void ParameterRemove(const char* parnam,element* el)
{
  const int ei=name_list_pos(parnam,el->def->par_names);
  if(ei > -1)
  {
    el->def->par_names->inform[ei]=0; // Turn off  by setting inform to 0   -- effective in save,  but still used in twiss; so really delete expression and turn value off
    const int ei_base=name_list_pos(parnam,el->def->par_names);
    double default_value=0;
    if(ei_base) default_value=el->base_type->def->par->parameters[ei_base]->double_value;  // el_par_value(parnam,el->base_type) cannot be used here, base element length may be zero
    command_parameter* cp=el->def->par->parameters[ei];
    if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " in " << el-> name << " parameter" << setw(12) << parnam
      << " value=" << setw(6) << cp->double_value << " set to default=" << setw(6) << default_value
      << " for " << setw(12) << parnam << " cp->expr=" << cp->expr << " and set expression to NULL" << EOL;
    cp->type = k_double;
    cp->double_value = default_value;
    cp->expr=NULL; // remove expression, -- without freeing space --  better memory leak than crash
    // if(cp->expr != NULL ) delete_expression(cp->expr); // delete expression --  resulted in crash in thick_bends/bend.madx  -- madX problem with poorly defined "objects"
  }
}

static string Check_command_parameter_consistence(const command* cmd)
{
  string EmptyStr="";
  if(cmd==NULL) return EmptyStr;
  command_parameter_list* cl=cmd->par;
  name_list* nl=cmd->par_names;
  if(cl==NULL || nl==NULL) return EmptyStr;
  if(cl->curr<1 || nl->curr<1) return EmptyStr;
  //-- at this point cmdpar and name_list exist and have at least one name,  see if they match
  vector<string> cl_names(cl->curr);
  for(int i=0;i < cl->curr; ++i) cl_names[i]=cl->parameters[i]->name;

  vector<string> nl_names(nl->curr);
  for (int i = 0; i < nl->curr; ++i) nl_names[i]=nl->names[i];
  unsigned int imax=cl_names.size();
  if(nl_names.size()>imax) imax=nl_names.size();
  int ierr=0;
  if(cl_names.size() != nl_names.size()) ierr=1;
  ostringstream ostr,ostr2;
  if(verbose_fl() || ierr) ostr2 << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " #cmdpar names=" << cl_names.size() << " #name_list names=" << nl_names.size() << EOL
	<< "    i         cmdpar name      name_list name" << EOL;
  for(unsigned int i=0; i<imax; ++i)
  {
	if(verbose_fl())
	{
	  ostr2 << setw(5) << i;
	  if(i<cl_names.size()) ostr2 << setw(20) << cl_names[i]; else ostr2 << "  NULL              ";
	  if(i<nl_names.size()) ostr2 << setw(20) << nl_names[i]; else ostr2 << "  NULL              ";
	  if(i<cl_names.size() && i<nl_names.size() && cl_names[i] !=nl_names[i]) ostr2 << " <----- not equal";
	  ostr2 << EOL;
	}
	if(i<cl_names.size() && i<nl_names.size() && cl_names[i] !=nl_names[i]) ierr++;
  }
  if(ierr)
  {
	ostr << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << "  *** ERROR *** " << " command " << cmd->name << " has inconsistent parameter_list and name_list names  ierr=" << ierr << EOL << ostr2.str() << EOL;
  }
  else if(verbose_fl()) ostr << "command " << cmd->name << " has consistent names" << EOL; // in this case no need to show ostr2
  // if(ierr) exit(ierr); // - exit in case of inconsistency - maybe useful for debugging after changes
  return ostr.str();
}

static string my_dump_command(const command* cmd)
{
  ostringstream ostr;
  if(cmd==NULL) ostr << " is NULL";
  else
  { // command defined in mad_elem.h, c-structures based on pointing to pointers, contains name_list par_names and command_parameter_list par
	ostr << "command: "    ; if(cmd->name)     ostr << cmd->name;   else ostr << " NULL";
	ostr << "  module: "   ; if(cmd->module)   ostr << cmd->module; else ostr << " NULL";
	ostr << "  group: "    ; if(cmd->group)    ostr << cmd->group;  else ostr << " NULL";
	ostr << "  stamp= " << cmd->stamp << "  link_type= " << cmd->link_type << "  mad8_type= " << cmd->mad8_type;
	ostr << "  #par_names "; if(cmd->par_names->curr)   ostr << cmd->par_names->curr; else ostr << " NULL";
	ostr << "  #par= "     ; if(cmd->par->curr)   ostr << cmd->group; else ostr << " NULL";
	ostr << EOL;
	ostr << "within command par_names:"; if(cmd->par_names) ostr << EOL << my_dump_name_list(cmd->par_names);        else ostr << " NULL" << EOL;
	ostr << "within command par:";       if(cmd->par)       ostr << EOL << my_dump_command_parameter_list(cmd->par); else ostr << " NULL" << EOL;
  }
  ostr << EOL;
  ostr << Check_command_parameter_consistence(cmd);
  return ostr.str();
}

static bool copy_cmd_par(const char* from_name,const char* to_name,const element* from_el,command* cmd) // copy parameter (expression) from an element to another command used in a new element, possible to change name
{
  command_parameter* to_param = clone_command_parameter( return_param(from_name,from_el) ); // clone to allow for changes in copy, uses my const clone_command_parameter which checks for NULL
  bool done=true,found=false;
  double value      =my_get_int_or_double_value(from_el           ,from_name,found);
  double default_val=my_get_int_or_double_value(from_el->base_type,from_name,found);
  const double eps=1.e-15; // used to check if strength is compatible with zero
  // if(to_param || found ) // use this to take anyway, even if on default value,  just for test, may result in  memory access outside program range
  if(to_param || (found && fabs(value-default_val)>eps) ) // expression defined or non-trivial value
  {
    if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " from_name=" << from_name << " to_name=" << to_name << " cloned to_param=" << to_param << " before strcpy " << my_dump_command_parameter(to_param);
    strcpy(to_param->name,to_name); // put to_name  to the cloned   command_parameter*
    if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " from_name=" << from_name << " to_name=" << to_name << " cloned to_param=" << to_param << " after  strcpy " << my_dump_command_parameter(to_param);
    add_cmd_parameter_clone(cmd,to_param,(to_name),1);
  }
  else
  {
    if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << from_name << " parameter not defined or on default values, do not copy " << my_dump_command_parameter(to_param) << '\n';
    done=false;
  }
  return done;
}

static void ParametersActiveOn(element* el) // turn active parameters on so that these values are written by save, where active means the expression exists or value different from default value,
{
  if(el && el->def && el->def->par)
  {
	command_parameter_list* cmdpar=el->def->par;
	for (int i = 0; i < cmdpar->curr; ++i)
	{
	  command_parameter* cmdi=cmdpar->parameters[i];
	  char* parnam=cmdi->name;
	  if( cmdi->expr &&  (strcmp(cmdi->expr->string, none) != 0) ) // turn on when expression defined and not none
	  {
		if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " in " << el-> name << " has expression, turn on " << parnam  << EOL;
		ParameterTurnOn(parnam,el);
	  }
	  else // check for non-default value
	  {
		const double eps=1.e-15; // used to check if strength is compatible with zero
		bool found=false;
		double default_val=my_get_int_or_double_value(el->base_type,parnam,found);
		if(found && fabs(cmdi->double_value-default_val)>eps)
		{
		  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " in " << el-> name
			<< " value=" << setw(6) << cmdi->double_value << " default=" << setw(6) << default_val << " has non trivial value, turn on " << parnam << EOL;
		  ParameterTurnOn(parnam,el); // turn this parameter on
		}
	  }
	}
  }
}

static string my_dump_element(const element* el)
{
  ostringstream ostr;
  ostr << setprecision(15) << EOL << "my_dump_element";
  if(el==NULL) ostr << " is NULL";
  else
  { // element defined in mad_cmd.h, c-structures based on pointing to pointers
	ostr << " name=" << el->name;
	if(el->base_type && el->base_type->name) ostr << " base_type=" << el->base_type->name;
	ostr << " stamp=" << el->stamp << " length=" << el->length << " parent name=" << el->parent->name;
	ostr << " def_type=" << el->def_type;
	if(el->def_type) ostr << " which means defined separately"; else ostr << " which means inside sequence";
	ostr << EOL;
	ostr << "within element " << my_dump_command(el->def);
  }
  return ostr.str();
}

static void Remove_All_Fringe_Field_Parameters(element* el)
{
#if __cplusplus >= 201103L // C++11 on
  vector<string> FringePar{"e1","e2","fint","fintx","h1","h2","hgap"};
#else
  vector<string> FringePar(7);
  FringePar[0]="e1";
  FringePar[1]="e2";
  FringePar[2]="fint";
  FringePar[3]="fintx";
  FringePar[4]="h1";
  FringePar[5]="h2";
  FringePar[6]="hgap";
  const char* FringePar_char[]={"e1","e2","fint","fintx","h1","h2","hgap"};
#endif
  const char* eltype=el->base_type->name;
  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " el name=" << el->name << " type" << eltype << " before remove : " << my_dump_element(el) << EOL;
  for(unsigned int i=0; i<FringePar.size(); ++i) ParameterRemove(FringePar[i].c_str(),el);
  if(kill_fringe_fl)
  {
    SetParameterValue("kill_ent_fringe",el,true,k_logical);
    SetParameterValue("kill_exi_fringe",el,true,k_logical);
    ParameterTurnOn("kill_ent_fringe",el); // turn writing on
    ParameterTurnOn("kill_exi_fringe",el); // turn writing on
  }
  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " el name=" << el->name << " type" << eltype << " after  remove : " << my_dump_element(el) << EOL;
}

static string my_dump_node(const node* node)
{
  ostringstream ostr;
  ostr << setprecision(15) << EOL << "my_dump_node";
  if(node==NULL) ostr << " is NULL";
  else
  {
	ostr << setprecision(15) << " node:";
	char pname[NAME_L] = "NULL", nname[NAME_L] = "NULL", from_name[NAME_L] = "NULL";
	if (node->previous != NULL) strcpy(pname, node->previous->name);
	if (node->next != NULL) strcpy(nname, node->next->name);
	if (node->from_name != NULL) strcpy(from_name, node->from_name);
	ostr << " name=" << left << setw(20) << node->name << right
	<< " occ=" << node->occ_cnt
	<< " node->base_name=" << left << setw(15) << node->base_name << right
	<< " from_name=" << left << setw(10) << from_name << right
	<< " at_value=" << setw(10) << node->at_value
	<< " position=" << setw(10) << node->position
	<< " previous=" << left << setw(15) << pname
	<< " next=" << setw(15) << nname << right
	<< " at_expr: ";
	if(node->at_expr) ostr << my_dump_expression(node->at_expr); else ostr << "NULL ";
	if(node->p_elem) ostr << my_dump_element(node->p_elem);
	if(node->cl!=NULL) for(int i=0; i< node->cl->curr; ++i) dump_constraint(node->cl->constraints[i]);
  }
  ostr << EOL;
  return ostr.str();
}

static string my_dump_sequence(const sequence* c_sequ,const int level)
{ // level 1 little info, 3 dump also nodes, 4 dump also elements
  ostringstream ostr;
  if(c_sequ==NULL) ostr << "sequence is NULL";
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
	ostr << " share=" << c_sequ->share << " nested=" << c_sequ->nested << " con_cnt=" << c_sequ->con_cnt << " stamp=" << c_sequ->stamp << " line=" << c_sequ->line << " add_pass=" << c_sequ->add_pass << " length=" << c_sequ->length << EOL;
	c_node = c_sequ->start;
    ostr << setprecision(15);
    double lastvalue=0;
	while(c_node != NULL)
	{
	  suml += c_node->length;
	  if (level > 2)
	  {
		ostr << my_dump_node(c_node);
		if (level > 3 && c_node->p_elem != NULL) ostr << my_dump_element(c_node->p_elem);
	  }
	  else if (level > 0 && strcmp(c_node->base_name, "drift") != 0)
	  {
		ostr << left << setw(20) << c_node->name << right
		<< " at_value=" <<  setw(8) << c_node->at_value
		<< " position=" <<  setw(8) << c_node->position
		<< " length="   << setw(17) << c_node->length;
		if(c_node->from_name) ostr << " from " << c_node->from_name;
		if(c_node->at_expr) ostr << " at_expr " << my_dump_expression(c_node->at_expr);
        if(c_node->at_expr) { double currentvalue=expression_value(c_node->at_expr, 2); ostr << " diff=" << currentvalue-lastvalue; lastvalue=currentvalue; }
		ostr << EOL;
	  }
	  if (c_node == c_sequ->end)  break;
	  c_node = c_node->next;
	}
	ostr << "===== sum of node length=" << setw(8) << suml << EOL;
	ostr << EOL;
  }
  return ostr.str();
}

static string my_get_cmd_expr_str(const command_parameter* cmd) // return the expression as string, if there is only a value, return the value as string
{
  string result="";
  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << my_dump_command_parameter(cmd);
  if(cmd)
  {
	if(cmd->expr  && cmd->expr->string) result=cmd->expr->string; // the expression is define as string, use this
	else // look for a value
	{
	  const double eps=1.e-15; // used to check if strength is compatible with zero
	  if( fabs(cmd->double_value)>eps ) // value defined and non-zero
	  {
		ostringstream ostr;
		if(cmd->double_value<0) ostr << "("; // enclose negative value in brackets
		ostr << cmd->double_value;
		if(cmd->double_value<0) ostr << ")"; // enclose negative value in brackets
		result=ostr.str();
	  }
	}
  }
  if(verbose_fl()) cout << " my_get_cmd_expr_str result=" << result << EOL;
  return result;
}

static expression* my_get_param_expression(const element* el,const char* parnam)// get a new copy of the expression for a parameter from an element, use the value as new expression if the expression was NULL
{
  const int ipar = name_list_pos(parnam,el->def->par_names);
  const command_parameter* cmdpar = NULL;
  if(ipar > -1)	cmdpar = el->def->par->parameters[ipar]; else return NULL; // pointer to the original length parameter
  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " for element " << setw(20) << el->name << " parameter " << setw(20) << parnam << " ipar=" << ipar << " my_dump_expression(cmdpar->expr):" << my_dump_expression(cmdpar->expr) << " cmdpar->double_value=" << cmdpar->double_value << EOL;
  command_parameter* cmdpar_copy = clone_command_parameter( cmdpar );  // copy of the original length parameter that can be modified
  // if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << "clone_command_parameter done" << EOL;
  if(cmdpar_copy->expr==NULL)
  { // use the value as new expression if the expression was NULL
	ostringstream ostr;
	ostr << setprecision(15) << cmdpar->double_value; // use the value as string
	cmdpar_copy->expr = new_expression(CopToNewC_String(ostr.str()),deco); // where deco is a global.  // cmdpar_copy->expr->value = cmdpar->double_value; // seems to have no effect and this not needed
	if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " cmdpar_copy->expr was NULL, create new expression from string " << ostr.str() << " now " << my_dump_expression(cmdpar_copy->expr) << EOL;
  }
  // if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << "done" << EOL;
  return cmdpar_copy->expr;
}

static bool thick_fl(const element* el) // true if the element has a thick parameter and if the value is positive, 0 otherwise
{
  const int thick_pos = name_list_pos("thick",el->def->par_names);
  return (thick_pos > -1 && el->def->par->parameters[thick_pos]->double_value>0);
}

static void dump_slices() // Loops over all current elements and prints the number of slices. Used for debug and info, uses global element_list
{
  unsigned int verbose=0;
  if(verbose_fl()) verbose=2;
  // verbose=3; // special for code development, shows all elements in element_list
  printf("++++++ dump_slices");
  if(verbose>1) printf("   verbose on, all elements are listed\n"); else printf("   only elements with non default selection (other than 1 thin) are shown\n");
  if(verbose>2) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " verbose=" << verbose << " list all elements" << EOL;
  printf("            name  #slices      derived from  #slices\n");
  int n_elem_with_slice=0,n_elem_with_slice_gt_1=0;
  for(int i=0; i< element_list->curr; ++i) // loop over element_list
  {
	element* el_i = element_list->elem[i];
	int el_i_slice_pos = name_list_pos("slice",el_i->def->par_names);
	if(el_i_slice_pos > -1) // element with slice number defined
	{
	  n_elem_with_slice++;
	  const int slices=el_i->def->par->parameters[el_i_slice_pos]->double_value;
	  int slices_parent=0;
	  char* parent_name=const_cast<char*>("no parent");
	  if(el_i->parent!=NULL) // look also at parent if existing
	  {
		slices_parent=el_i->parent->def->par->parameters[el_i_slice_pos]->double_value;
		parent_name=el_i->parent->name;
	  }
	  if(slices>1) n_elem_with_slice_gt_1++;
	  if(verbose_fl() || slices !=1 || thick_fl(el_i)) // print all with verbose. with debug skip elements with default selection thin 1 slice
	  {
		printf(" %15s %2d",el_i->name,slices);
		if(thick_fl(el_i)) printf(" thick"); else printf(" thin ");
		if(el_i != el_i->parent)
		{
		  printf("%18s %2d",parent_name,slices_parent); // show also parent if not same as child
		  if(thick_fl(el_i->parent)) printf(" thick"); else printf(" thin ");
		}
		printf("\n");
	  }
	}
	else if(verbose>2) cout << setw(16) << el_i->name << EOL; // no slice number defined, just give the name
  }
  if(verbose>1) cout << "       general option thin_foc=" << thin_foc_fl() << EOL; // global option like debug or verbose, not element specific, still print here for info
  printf("------ end of dump slices. There were %4d elements, %3d with slice numbers and %2d with slice numbers>1\n\n",element_list->curr,n_elem_with_slice,n_elem_with_slice_gt_1);
}

static void force_consistent_slices() // hbu 10/2005 loop over all elements and check that #slices of child and parent agree,  if not, use the maximum for both
{
  for(int i=0; i< element_list->curr; ++i) // loop over element_list
  {
	element* el_i = element_list->elem[i];
	int el_i_slice_pos = name_list_pos("slice",el_i->def->par_names);
	if(el_i_slice_pos > -1 && el_i->parent!=NULL && el_i != el_i->parent )
	{
	  command_parameter*  child=el_i->def->par->parameters[el_i_slice_pos];
	  command_parameter* parent=el_i->parent->def->par->parameters[el_i_slice_pos];
	  int slices=child->double_value;
	  int slices_parent=parent->double_value;
	  if(slices != slices_parent)
	  {
		if(slices>slices_parent) slices_parent=slices; else slices=slices_parent;
		child->double_value=slices;
		parent->double_value=slices_parent;
		int el_i_thick_pos = name_list_pos("thick",el_i->def->par_names);
		if(el_i_thick_pos > -1) el_i->parent->def->par->parameters[el_i_thick_pos]->double_value = el_i->def->par->parameters[el_i_thick_pos]->double_value; // copy thick flag from child to parent
	  }
	}
  }
  if (debug_fl()) { printf("end of force_consistent_slices\n"); dump_slices(); }
}

static int get_slices_from_elem(const element* elem)
{
  int elem_slice_pos=0,slices=1;
  if( (elem_slice_pos = name_list_pos("slice",elem->def->par_names)) > -1 ) slices=elem->def->par->parameters[elem_slice_pos]->double_value;
  return slices;
}

static char* make_thin_name(const char* e_name,const int slice) // make a node name from element name and slice number
{ // example     e_name=mqxa.1r1 slice=1 result=mqxa.1r1..1
  ostringstream ostr;
  ostr << e_name << ".." << slice;
  return CopToNewC_String(ostr.str());
}

static command_parameter*
scale_and_slice(command_parameter* kn_param,const command_parameter* length_param,const int slices,const int angle_conversion,const int kl_flag) // multiply the k by length and divide by slice
{
  int last_non_zero=-1;
  if (kn_param == NULL) return NULL;
  const double eps=1.e-15; // used to check if strength is compatible with zero

  for (int i=0; i<kn_param->expr_list->curr; ++i)
  {
	expression* kn_i_expr = kn_param->expr_list->list[i];
	double kn_i_val  = kn_param->double_array->a[i];
	if ((kn_i_expr != NULL && zero_string(kn_i_expr->string)==0) || fabs(kn_i_val)>eps )
	{
	  last_non_zero=i;
	  if (kl_flag == 0 && (angle_conversion==0||i>0)) //hbu apply the angle_conversion==0 check only to zero order multipole
	  {
		if ((length_param->expr) || (kn_i_expr)) kn_i_expr = compound_expr(kn_i_expr,kn_i_val,const_cast<char*>("*"),length_param->expr,length_param->double_value); // multiply expression with length
		else kn_i_val *= length_param->double_value; // multiply value with length
	  }
	  if (slices > 1) // give the correct weight by slice (multiply with the inverse of the number of slices)
	  {
		if (kn_i_expr) kn_i_expr = compound_expr(kn_i_expr,kn_i_val,const_cast<char*>("*"),NULL,1./slices);
		else kn_i_val *= 1./slices;
	  }
	  if(verbose_fl()) { printf("verbose %s %s line %d  kn_i_val=%f  kl_flag=%d\n",__FILE__,__FUNCTION__,__LINE__,kn_i_val,kl_flag); if(kn_i_expr) cout << my_dump_expression(kn_i_expr) << EOL; }
	}
	if(kn_i_expr) kn_param->expr_list->list[i] = kn_i_expr;
	kn_param->double_array->a[i] = kn_i_val;
  } // for i ..
  if (last_non_zero==-1)
  {
	delete_command_parameter(kn_param);
	kn_param=NULL;
  }
  return kn_param;
}

static void add_lrad(command* cmd,const command_parameter* length_param,const int slices)
{
  command_parameter* l_par;
  if(length_param)
  {
	add_cmd_parameter_new(cmd,0.,"l",1); // new parameter l with value of 0
	l_par = cmd->par->parameters[cmd->par->curr] = clone_command_parameter(length_param); // keep what was l
	strcpy(l_par->name,"lrad"); // but rename to lrad and slice :
	if (slices > 1) // divide numbers or expressions by the number of slices
	{
	  if (l_par->expr) l_par->expr = compound_expr(l_par->expr,0.,const_cast<char*>("/"),NULL,slices);
	  else l_par->double_value /= slices;
	}
	add_to_name_list(const_cast<char*>("lrad"),1,cmd->par_names);
	cmd->par->curr++;
  }
}

static node* new_marker(const node* thick_node,const double at,expression *at_expr) // create a new marker element and node with position given by at_expr
{
  
#if __cplusplus >= 201103L // C++11 on
  vector<const char*> ParList{"at","kmax","kmin","polarity","calib","type","apertype","aperture","aper_tol","mech_sep","mech_sep","v_pos","from"};
#else // old C++ compiler
  vector<const char*> ParList(12);
  ParList[0]="at";
  ParList[1]="kmax";
  ParList[2]="kmin";
  ParList[3]="polarity";
  ParList[4]="calib";
  ParList[5]="type";
  ParList[6]="apertype";
  ParList[7]="aperture";
  ParList[8]="aper_tol";
  ParList[9]="mech_sep";
  ParList[10]="v_pos";
  ParList[11]="from";
#endif
  node* node=NULL;
        element* elp=thick_node->p_elem;
  const element* thick_el=thick_node->p_elem;
  if (elp)
  {
	int imarker_pos = name_list_pos("marker", defined_commands->list);
	command* p = defined_commands->commands[imarker_pos];
	const int mx=25; // maximum number of par names and parvalues
	command* cmd = new_command(p->name, mx, mx, p->module, p->group, p->link_type,p->mad8_type);
    for(unsigned int i=0; i<ParList.size(); ++i) add_cmd_parameter_clone(cmd,return_param_recurse(ParList[i], thick_el),ParList[i],1);
	element* elem = make_element(elp->name,const_cast<char*>("marker"), cmd,-1);
	if ( verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " new marker element " << my_dump_element(elem) << EOL;
	node = new_elem_node(elem, thick_node->occ_cnt);
	strcpy(node->name, thick_node->name);
	node->occ_cnt = thick_node->occ_cnt;
	node->at_value = at;
	if (at_expr) node->at_expr = clone_expression(at_expr);
	node->from_name = thick_node->from_name;
  }
  else
  {
	fatal_error("This is not an element ",thick_node->name);
  }
  return node;
}

static expression* curved_from_straight_length(const element* rbend_el)
{
  expression* l_rbend_expr = my_get_param_expression(rbend_el,const_cast<char*>("l")); // get expression or create new from constant
  expression* l_sbend_expr = NULL;
  if(rbarc_fl() && l_rbend_expr ) // increase the straight rbend length to sbend length
  { // Mad-X very confusing for RBEND, "l" parameter   el_par_value("l","rbend") is the shorter straight length,   val = l * angle / (two * sin(angle/two));     with rbarc on as done by default
    // this is also shown in twiss      node and element length give always the curved length
    // in going from rbend to sbend, this correction must be applied   if the "l" expression is used,    not for the value
    const string anglestr    =my_get_cmd_expr_str( return_param_recurse(const_cast<char*>("angle"),rbend_el) );
    const string rat = "("+anglestr+")*0.5/sin(("+anglestr+")/2)"; // L_sbend / L_rbend
    expression* rat_expr = new_expression(CopToNewC_String(rat),deco);
    // try status=0 or 1 to update
    l_sbend_expr = compound_expr(l_rbend_expr,0,const_cast<char*>("*"),rat_expr,0); // this also updates the value
    if(verbose_fl())
    {
      bool found=false;
      double straight_length=my_get_int_or_double_value(rbend_el,"l",found);
      cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__  << " " << rbend_el->name << " rbarc on, increase rbend straight length expression of value " << straight_length << " to curved sbend length  using anglestr=" << anglestr
      << " updated l_sbend_expr " << my_dump_expression(l_sbend_expr) << " value should now be same as the curved rbend_el->length=" << rbend_el->length << EOL;
    }
  }
  else l_sbend_expr=l_rbend_expr;
  return l_sbend_expr;
}

static element* sbend_from_rbend(const element* rbend_el,const bool MakeDipedge)
{
  // go from rbend to sbend
  // just changing the base_name did not work - seems also to change all parents, even with clone
  // if done on parent, then the next child will think it is already sbend and stop conversion
  // so rather go for a completely new sbend and copy only what is needed. Do not copy the parent rbend
  // new_element in older makethin did not this as needed here, rather construct a new command and then use make_element similar to what is done here for other new elements like dipege
  // leave the rbend thick_elem as they were, just do not place them
  unsigned int verbose=0;
  if(debug_fl()) verbose=1;
  if(verbose_fl()) verbose=2;

  const string sbend_name=string(rbend_el->name)+"_sbend";

  command* sbend_def = defined_commands->commands[ name_list_pos("sbend", defined_commands->list) ];
  if(verbose>1)
  {
	cout << EOL << EOL;
	cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__;
	cout << " for rbend_el " << rbend_el->name;
	if(verbose>2) cout << my_dump_element(rbend_el);
	cout << " sbend_name=" << sbend_name
	<< " def->name=" << sbend_def->name
	<< " def->module=" << sbend_def->module  // element
	<< " sbend_def->par->curr=" << sbend_def->par->curr
	<< EOL;
	if(verbose>2) cout << " where sbend defined :" << my_dump_command(sbend_def) << EOL;   //  count the number of parameters,  seems  47 ?
  }

  command* sbend_cmd=NULL;
  const bool CloneCommand=false; //---  with false  no diff thin and thin2
  if(CloneCommand)
  {
    sbend_cmd = clone_command(rbend_el->def); // start from clone of rbend defining command
    if(verbose>1) cout  << EOL << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << "  CloneCommand=" << CloneCommand << " " << my_dump_command(sbend_cmd) << EOL;
  }
  else // also start defining command from scratch -  works, except that difficult to get extra parameters  thick, kill_ent_fringe, kill_exi_fringe turned on
  {
    const int mx=sbend_def->par->curr; // maximum number of par names and parvalues -- take from definition,   getting mx=47  1/11/2014
    sbend_cmd = new_command(CopToNewC_String(string(sbend_def->name)+"_cmd"), mx, mx, sbend_def->module, sbend_def->group, sbend_def->link_type,sbend_def->mad8_type); // new command, used here to define sbend
    if(verbose>1) cout  << EOL << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " mx=" << mx << " new sbend_cmd " << my_dump_command(sbend_cmd) << EOL;
  }

#if __cplusplus >= 201103L // C++11 on
  vector<string> LogParListToCopy{"thick","kill_ent_fringe","kill_exi_fringe"};
#else // old C++ compiler
  vector<string> LogParListToCopy(3);
  LogParListToCopy[0]="thick";
  LogParListToCopy[1]="kill_ent_fringe";
  LogParListToCopy[2]="kill_exi_fringe";
#endif

  if(rbend_el && rbend_el->def && rbend_el->def->par)
  {
	command_parameter_list* cmdpar=rbend_el->def->par;
	for (int i = 0; i < cmdpar->curr; ++i)
	{
	  int inform=0; // for  add_cmd_parameter_new/add_to_name_list
	  command_parameter* cmdi=cmdpar->parameters[i];
      double value=cmdi->double_value;
	  char* parnam=cmdi->name;
	  if( cmdi->expr && (strcmp(cmdi->expr->string, none) != 0) ) // turn on when expression defined and not none
	  {
		if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " in " << rbend_el-> name << " has expression, use this " << parnam  << EOL;
		add_cmd_parameter_clone(sbend_cmd, return_param_recurse(parnam,rbend_el),parnam,1);
	  }
      else if(NameIsInList(string(parnam),LogParListToCopy)) // in LogParListToCopy   add to list
      {
        add_cmd_parameter_new(sbend_cmd,value,parnam,inform); // by default double, use SetParameterValue with k_logical for bool
      }
	  else // copy value,   just copying the non-default values enough for save, twiss instead seems to also require default values
	  {
		const double eps=1.e-15; // used to check if strength is compatible with zero
		bool found=false;
		double default_val=my_get_int_or_double_value(rbend_el->base_type,parnam,found); // return the default values, works for int and double parameters
        if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " rbend " << rbend_el->name << " parnam " << parnam << " " << sbend_name
          << " cmdi=" << cmdi << " value=" << value << EOL;
        if( string(parnam) == string("l") && rbarc_fl())
        {
          value = el_par_value(parnam,rbend_el); // special case rbend with option("rbarc"), get increased arc length value from el_par_value
          if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " rbarc on, use increased length value=" << value << EOL;
        }
		if(found) //  && fabs(value-default_val)>eps ) // adding only non-default values ok for save, but not for twiss, make all parameters, even if default value
		{
		  if (fabs(value-default_val)>eps ) inform=1; // different from default, mark to be written in safe
		  if(verbose_fl())
		  {
			cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " in " << rbend_el-> name
			<< " value=" << setw(6) << value << " default=" << setw(6) << default_val << " use this " << parnam;
			if(inform) cout << " differs from default, set inform=" << inform;
			cout << EOL;
		  }
		  add_cmd_parameter_new(sbend_cmd,value,parnam,inform);
		}
	  }
	}
  }
  // --   at this point, the rbend length should already be in sbend_cmd,   needs only modification in special case
  if( verbose>1 ) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << my_dump_command(sbend_cmd);
  if(rbarc_fl())
  {
	expression* l_sbend_expr=curved_from_straight_length(rbend_el); // use this modified length expression in sbend_cmd
	int il=name_list_pos(const_cast<char*>("l"),sbend_cmd->par_names); // parameter 0
	if(il > -1) sbend_cmd->par->parameters[il]->expr=l_sbend_expr;
	if( verbose>1 ) cout << EOL << EOL << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " after increase of rbend length now l_sbend_expr : "
      << my_dump_expression(l_sbend_expr)
      << " sbend_cmd : " << my_dump_command(sbend_cmd);
  }
  if(verbose>1) cout << EOL << EOL << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << my_dump_command(sbend_cmd) << EOL;
  if(verbose>1) cout << EOL << EOL << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " just before element *sbend_el=make_element  sbend_name=" << sbend_name << EOL;
  element *sbend_el=make_element(const_cast<char*>(sbend_name.c_str()),const_cast<char*>("sbend"),sbend_cmd,2);

  for(unsigned int i=0; i<LogParListToCopy.size(); ++i)
  {
    const char* parnam=LogParListToCopy[i].c_str();
    SetParameterValue(parnam,sbend_el,true,k_logical);
    if(MakeDipedge) ParameterTurnOn(parnam,sbend_el); //-- so that thick=true is written  to signal this for thick tracking
  }

  if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__  << " now sbend_el=" << sbend_el << " " << my_dump_element(sbend_el) << " compare this to the original rbend " << my_dump_element(rbend_el) << EOL;
  return sbend_el;
} // sbend_from_rbend

static void add_node_at_end_of_sequence(node* node,sequence* sequ) // position in thin sequence defined with at_value and from_name
{
  if (sequ->start == NULL) // first node in new sequence
  {
	sequ->start = node;
	node->next = NULL;
	node->previous = NULL;
  }
  else // add to end
  {
	sequ->end->next = node;
	node->previous  = sequ->end;
  }
  sequ->end = node;

  // node->at_expr=NULL; // use this to only write at values, no expressions

  if(verbose_fl())
  {
	cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << left << setw(20) << node->name << " " << setw(20) << node->base_name << right
	<< " position=" << setw(10) << node->position  << " at_value=" << setw(10) << node->at_value;
	if(node->at_expr)   cout << " " << my_dump_expression(node->at_expr);
	if(node->from_name) cout << " from " << setw(5) << node->from_name; else cout << "           ";
	cout << " length="   << setw(10) << node->length  << " in " << sequ->name << EOL;
  }
  add_to_node_list(node, 0, sequ->nodes);
  return;
}

static void add_half_angle_to(const element* rbend_el,element* to_el,const char* to_parm) // get half surface angle of rbend, add to e1 or e2 of dipedge or sbend
{
  if(rbend_el && to_el)
  {
    expression* half_angle_expr  = scale_expr(my_get_param_expression(rbend_el,"angle"),0.5); // angle*0.5
    command_parameter* to_param = return_param_recurse(const_cast<char*>(to_parm),to_el); // get param from element, may not exist, use here the non const version of return_param_recurse, to modify to_el
    if(to_param) // modify the existing parameter in to_el
    {
      if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " original to_param " << my_dump_command_parameter(to_param) << EOL;
      to_param->expr = compound_expr(to_param->expr,0,const_cast<char*>("+"),half_angle_expr,0);
      if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << "    now  to_param " << my_dump_command_parameter(to_param) << EOL;
    }
    else // param not yet in to_el, start from parameter definition
    {
      int ipar = name_list_pos(to_parm,to_el->def->par_names);
      if(ipar > -1) // already in name_list
      {
        to_param = clone_command_parameter( to_el->def->par->parameters[ipar] );  // copy of the original length parameter that can be modified
        to_el->def->par->parameters[ipar]->expr=half_angle_expr;
        ParameterTurnOn(to_parm,to_el);
        if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " use existing to_param from ipar= " << ipar
          << " to_param=" << to_param
          << " to_el->def->par->parameters[ipar]->expr=" << to_el->def->par->parameters[ipar]->expr << EOL;
      }
      else // not in name_list_pos
      {
        ipar = name_list_pos(to_parm,to_el->base_type->def->par_names); // parameter in the definition, must always exist
        if(ipar > -1)
        {
          to_param = clone_command_parameter( to_el->base_type->def->par->parameters[ipar] );  // copy of the original length parameter that can be modified
          if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " for element " << setw(20) << to_el->name << " parameter " << setw(20) << to_parm << " my_dump_expression(to_param->expr):" << my_dump_expression(to_param->expr) << " to_param->double_value=" << to_param->double_value << EOL;
        }
        else
        {
          if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " *** error ***<  element " << setw(20) << to_el->name << " of type " << to_el->base_type->name << " has no parameter" << to_parm << EOL;
          return;
        }
        to_param->expr=half_angle_expr; // set expression, which is NULL from definition, to half_angle_expr
        if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << "    now  to_param " << my_dump_command_parameter(to_param) << EOL;
        add_cmd_parameter_clone(to_el->def,to_param,to_parm,1); // add new parameter to element, increases the name_list
      }
    }
  }
}

static element* create_bend_dipedge_element(element* thick_elem,const bool Entry) // using Dipedge  http://mad.web.cern.ch/mad/madx.old/Introduction/dipedge.html, also return sbend_el in case of rbend
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
  unsigned int verbose=0;
  if(debug_fl()) verbose=1;
  if(verbose_fl()) verbose=2;
  const bool IsRbend= string(thick_elem->base_type->name) == string("rbend");

  vector<string> CheckParams; //-- parameters which can directly be copied from bend to dipedge if value different from default. Instead e1, e2, fint, fintx, h1, h2 need extra work
  CheckParams.push_back("polarity");
  CheckParams.push_back("tilt");
  //old CheckParams.push_back("h");
  CheckParams.push_back("hgap");
  CheckParams.push_back("mech_sep");
  CheckParams.push_back("v_pos");
  CheckParams.push_back("magnet");
  CheckParams.push_back("model");
  CheckParams.push_back("method");
  CheckParams.push_back("exact");
  CheckParams.push_back("nst");

  // verbose=3; // extra debug - print element definitions

  element* dipedge=NULL;
  if (thick_elem)
  {
	string dipedge_name=string(thick_elem->name);
	command* dipedge_def = defined_commands->commands[ name_list_pos("dipedge", defined_commands->list) ]; // access to dipedge structure as defined in mad_dict.h
	if(Entry) dipedge_name+="_den"; else dipedge_name+="_dex"; // dipedge entry or exit
	if(verbose>1)
	{
	  cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__;
	  if(Entry) cout << " at entry"; else cout << " at  exit";
	  cout << " for thick_elem " << thick_elem->name;
	  if(verbose>2) cout << my_dump_element(thick_elem) << EOL;
	  cout << " dipedge_name=" << dipedge_name
	  << " def->name=" << dipedge_def->name  // dipedge      mad8_type=33
	  << " def->module=" << dipedge_def->module  // element
	  << EOL;
	  if(verbose>2) cout << " where dipedge defined :" << my_dump_command(dipedge_def) << EOL;   //  count the number of parameters,  seems  47 ?
	}
	int mx=dipedge_def->par->curr; // maximum number of par names and parvalues -- take from definition
	if(verbose>1) cout << " dipedge_def->par->curr=" << dipedge_def->par->curr << EOL;
	string dipedge_cmd_name=string(dipedge_def->name);
	if(Entry) dipedge_cmd_name+="_l_"; else dipedge_cmd_name+="_r_";
	dipedge_cmd_name+="cmd";
	command* dipedge_cmd = new_command(CopToNewC_String(dipedge_cmd_name), mx, mx, dipedge_def->module, dipedge_def->group, dipedge_def->link_type,dipedge_def->mad8_type); // new command, used here to define dipedge
	const double eps=1.e-15;
	const bool generate_h=true;
	
	expression* l_par_expr=NULL;
	if(generate_h) // from length and angle
	{
	  char* parnam=const_cast<char*>("l");
	  const int ipar = name_list_pos(parnam,thick_elem->def->par_names);
	  command_parameter* cmdpar = NULL;
	  if(ipar > -1)
	  {
        cmdpar = thick_elem->def->par->parameters[ipar]; // pointer to the original length parameter
		if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " for element " << setw(20) << thick_elem->name << " parameter " << setw(20) << parnam << " my_dump_expression(cmdpar->expr):" << my_dump_expression(cmdpar->expr) << " cmdpar->double_value=" << cmdpar->double_value << EOL;
		command_parameter* cmdpar_copy = clone_command_parameter( cmdpar );  // copy of the original length parameter that can be modified
		l_par_expr=cmdpar_copy->expr;
	  }
	  if(l_par_expr) // length expression
	  {
		if(l_par_expr && IsRbend && rbarc_fl() ) l_par_expr=curved_from_straight_length(thick_elem); // increase straight length to curved length
	  }
	  else // only length value, use it and put in expression
	  {
		double l_par_value= el_par_value(const_cast<char*>("l"),thick_elem); // also does straight length to curved length conversion if needed
		if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__  << " no length expression, only value=" << l_par_value << EOL;
		int ipar = name_list_pos(const_cast<char*>("l"),thick_elem->def->par_names);
		if(ipar > -1)
		{
		  ostringstream ostr;
		  ostr << setprecision(15) << l_par_value; // use the value as string
		  l_par_expr = new_expression(CopToNewC_String(ostr.str()),deco); // where deco is a global.
		}
	  }
	  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " l_par_expr= " << my_dump_expression(l_par_expr) << EOL;
	  expression* angle_par_expr = my_get_param_expression(thick_elem,const_cast<char*>("angle"));
	  command_parameter* hparam=new_command_parameter(CopToNewC_String("h"),k_double);
	  hparam->expr=compound_expr(angle_par_expr,0.,const_cast<char*>("/"),l_par_expr,0); // this also updates the value
	  add_cmd_parameter_clone(dipedge_cmd,hparam,const_cast<char*>("h"),1);
	}

    if(dipedge_h1_h2_fl)
    {
      if(Entry) copy_cmd_par("h1","h1",thick_elem,dipedge_cmd); // at entry, copy h1 from thick bend as dipedge h1
      else      copy_cmd_par("h2","h2",thick_elem,dipedge_cmd); // at  exit, copy h2 from thick bend as dipedge h2
    }

    if(Entry) copy_cmd_par("e1","e1",thick_elem,dipedge_cmd); // at entry, copy e1 from thick sbend as dipedge e1
    else      copy_cmd_par("e2","e1",thick_elem,dipedge_cmd); // at  exit, copy e2 from thick sbend as dipedge e1
    
    if(Entry)
    {
      copy_cmd_par("fint","fint",thick_elem,dipedge_cmd); // at entry, copy fint from thick bend as dipedge fint
    }
    else // Exit
    {
      const command_parameter *fintxparam = return_param("fintx",(const element*)thick_elem); // use my const version of return_param to check the presence of fintx in thick_elem
      if(fintxparam!=NULL) copy_cmd_par("fintx","fint",thick_elem,dipedge_cmd); // use fintx if given, no check on value, allows for fintx=0 to turn off exit
      else copy_cmd_par("fint" ,"fint",thick_elem,dipedge_cmd); // use fint
      if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << my_dump_command_parameter(fintxparam) << EOL;
    }

    if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << my_dump_command(dipedge_cmd) << EOL;

    bool found=false;

    for(unsigned int i=0; i<CheckParams.size(); ++i) // copy other nontrivial parameters given in CheckParams from thick bend  -- only gets here with nontrivial e1 or e2     -- otherwise already returned NULL before
	{
	  const string ParNam=CheckParams[i];
	  char* parnam=CopToNewC_String(ParNam.c_str()); // use this,    char* parnam=const_cast<char*>(ParNam.c_str());   resulted in name_list corruption
	  command_parameter* this_param = return_param(parnam,thick_elem);
	  double value      =my_get_int_or_double_value(thick_elem           ,parnam,found);
	  double default_val=my_get_int_or_double_value(thick_elem->base_type,parnam,found);
	  if(this_param || (found && fabs(value-default_val)>eps) ) // expression defined or non-trivial value
	  {
		add_cmd_parameter_clone(dipedge_cmd,this_param,parnam,1); // use this parameter (expression or value) for dipedge
		if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << thick_elem->name << "        use parameter " << setw(12) << parnam << " for dipedge this_param=" << this_param << EOL;
	  }
	  else
	  {
		if(verbose>2) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << thick_elem->name << " do not use parameter " << setw(12) << parnam << " for dipedge this_param=" << this_param << EOL;
	  }
	}
	dipedge=make_element(const_cast<char*>(dipedge_name.c_str()),const_cast<char*>("dipedge"), dipedge_cmd,-1); // make the element and put it in global element_list, using the command dipedge_cmd, -1 means avoid warnings,  1 means delete and warn, 2 means warn and ignore if already present  see  add_to_el_list  mad_elem.c
    
    if(verbose>2)
	{
	  cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << setw(12) << thick_elem->name << " after ";
	  if(Entry) cout << "e1"; else cout << "e2";
	  cout << " remove" << EOL << my_dump_name_list(thick_elem->def->par_names) << EOL << my_dump_command_parameter_list(thick_elem->def->par) << EOL;
	}
  }
  if(verbose>2) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " at end of create_bend_dipedge_element " << my_dump_element(dipedge) << " dipedge address=" << dipedge << EOL;
  return dipedge;
}

static void create_and_place_bend_node(const node* thick_node,expression *at_expr,element* slice_elem, sequence* sequ,const double rel_shift,const bool IsRbend) // node to place the dipedge
{
  if(thick_node->p_elem)
  {
	node*  slice_node = new_elem_node(slice_elem, thick_node->occ_cnt);
	slice_node->from_name = buffer(thick_node->from_name);
	slice_node->occ_cnt = thick_node->occ_cnt;   //  occ_cnt is the occurance count, take from thick element
	slice_node->at_value = thick_node->at_value;
	slice_node->length = thick_node->length;
	slice_node->from_name = thick_node->from_name;
	
	if(slice_elem->base_type && slice_elem->base_type->name && (strcmp(slice_elem->base_type->name,"dipedge") == 0) ) slice_node->length=0; // set dipedge length = 0   (otherwise copied from thick bend)

	if(verbose_fl())
	{
	  cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " make node " << slice_elem->name << " " << slice_elem->base_type->name << " original thick " << thick_node->p_elem->name << my_dump_element(slice_elem) << EOL
	  << " slice_node->from_name=";
	  if(slice_node->from_name) cout << slice_node->from_name; else cout << "NULL ";
	  cout << " thick_node->at_value=" << thick_node->at_value;
	  cout << " slice_node->at_value=" << slice_node->at_value << EOL;
	}
	expression* length_param_expr=my_get_param_expression(thick_node->p_elem,const_cast<char*>("l")); // get expression or create new from constant
    if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " length_param_expr " << my_dump_expression(length_param_expr) << EOL;
	if(IsRbend) length_param_expr = curved_from_straight_length( thick_node->p_elem);
    if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " length_param_expr " << my_dump_expression(length_param_expr) << EOL;
	
	double at = thick_node->at_value;
	slice_node->at_expr = compound_expr(at_expr, at,const_cast<char*>("+"),  scale_expr(length_param_expr,rel_shift),  0 ); // this also updates the value
	
	if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__
	  << " length_param_expr " << my_dump_expression(length_param_expr)
	  << " at=" << at
	  << " slice_node->at_expr=" << my_dump_expression(slice_node->at_expr) << EOL;
	
	// slice_node->at_expr = at_expr;
	slice_node->at_value  = at;
	add_node_at_end_of_sequence(slice_node,sequ);	
  }
  else
  {
	fatal_error("This is not an element ",thick_node->name);
  }
}

static void place_thick_slice(element* thick_elem,const node* node, sequence* to_sequ, element* slice_elem, const int i, const string& thin_style) // make nodes for the _s, _b  pieces  and place them in the sequence
{
  const int n = get_slices_from_elem(thick_elem);
  const bool IsRbend= string(thick_elem->base_type->name) == string("rbend");
  SliceDistPos SP(n, thin_style==string("teapot") );
  if(verbose_fl())
  {
    cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " n=" << n << " start from thick_elem " << thick_elem->name;
    if(IsRbend) cout << " rbend -- careful thick_elem->length=" << thick_elem->length << " is the curved length, longer than the straight length in element";
    cout << my_dump_element(thick_elem); SP.Print();
  }

  expression* l_par_expr;
  if(IsRbend  && rbarc_fl() ) l_par_expr=curved_from_straight_length(thick_elem);
  else l_par_expr=my_get_param_expression(thick_elem,const_cast<char*>("l")); // with this l_par_expr should not be NULL

  expression* at_expr = clone_expression(node->at_expr);
  double at = node->at_value;

  double rel_shift;	
  if(i==0)       rel_shift=-0.5 + SP.delta/2.; // entry
  else if(i==n)  rel_shift= 0.5 - SP.delta/2.; // exit
  else           rel_shift=-0.5 + SP.delta + (i-0.5)*SP.Delta; // body

  at_expr = compound_expr(at_expr,at,const_cast<char*>("+"),  scale_expr(l_par_expr,rel_shift),  0 ); // this also updates the value
  struct node* thick_node = new_elem_node(slice_elem, node->occ_cnt);
  thick_node->from_name = buffer(node->from_name);
  thick_node->at_expr   = at_expr;
  thick_node->at_value  = at;
  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " place " << slice_elem->name << " using at_expr where " << my_dump_expression(at_expr) << " at_value=" << at << EOL;
  add_node_at_end_of_sequence(thick_node,to_sequ); // add the thick node to the sequences
  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " dump_node(thick_node)   :" << EOL << my_dump_node(thick_node) << " done with place_thick_slice " << i << EOL;
}


static sequence* slice_sequence(const string& thin_style,sequence* thick_sequ) // make recursively a sliced sequence out of the thick_seque
{
  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " thin_style=\"" << thin_style << "\"" << EOL;

  sequence* thin_sequ;
  static SequenceList sliced_seqlist;
  if ((thin_sequ=sliced_seqlist.get_sequ(thick_sequ))) return thin_sequ; // do nothing if the sequence was already sliced

  char name[128];
  strcpy(name,thick_sequ->name);
  fprintf(prt_file, "makethin: slicing sequence : %s\n",name);

  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " my_dump_sequence thick_sequ " << my_dump_sequence(thick_sequ,2) << EOL; // dump level 2, without nodes/elements

  thin_sequ = new_sequence(name, thick_sequ->ref_flag);
  thin_sequ->start = NULL;
  thin_sequ->share = thick_sequ->share;
  thin_sequ->nested = thick_sequ->nested;
  thin_sequ->length = sequence_length(thick_sequ);
  thin_sequ->refpos = buffer(thick_sequ->refpos);
  thin_sequ->beam = thick_sequ->beam;
  if (thin_sequ->cavities != NULL)  thin_sequ->cavities->curr = 0;
  else thin_sequ->cavities = new_el_list(100);
  if (thin_sequ->crabcavities != NULL)  thin_sequ->crabcavities->curr = 0;
  else thin_sequ->crabcavities = new_el_list(100);

  SeqElList theSeqElList(name,thin_style,thick_sequ,thin_sequ,thick_sequ->start);
  while(theSeqElList.thick_node != NULL) // loop over current sequence, the nodes are added to the sequence in add_node_at_end_of_sequence()
  {
	theSeqElList.slice_node();
	if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << my_dump_node( theSeqElList.thick_node );
	if (theSeqElList.thick_node == thick_sequ->end)  break;
	theSeqElList.thick_node = theSeqElList.thick_node->next;
  }
  thin_sequ->end->next = thin_sequ->start;

  if(verbose_fl())  cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << EOL << EOL << EOL << " my_dump_sequence thin_sequ " << my_dump_sequence(thin_sequ,2)  << EOL; // dump level 2, without nodes/elements

  int pos=0;
  if ((pos = name_list_pos(name, sequences->list)) < 0) // move the pointer in the sequences list to point to our thin sequence
  {
	fatal_error("unknown sequence sliced:", name);
  }
  else
  {
	sequences->sequs[pos]= thin_sequ; // pointer moved ok, delete_sequence(thick_sequ)
  }
  sliced_seqlist.put_sequ(thick_sequ); // Slicing done for this sequence. Add to list of sequences sliced
  if(debug_fl()) theSeqElList.Print(); // final list
  return thin_sequ;
} // slice_sequence

static int set_selected_elements() //  result in global  element_list     used in dump_slices, force_consistent_slices
{
  if (current_sequ == NULL || current_sequ->ex_start == NULL) // check that there is an active sequence, would otherwise crash in get_ex_range
  {
	warning("makethin selection without active sequence,", "ignored");
	return 1;
  }
  // select, flag=makethin, clear;   only resets the selection commands by  slice_select->curr=0  in mad_select.c
  // the element_list ist not cleared, reset it here, before evaluation of the selection for makethin
  for(int j=0; j< element_list->curr; ++j) // loop over element_list
  {
	element* el_j = element_list->elem[j];
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
	command_parameter_list* pl = slice_select->commands[i]->par;
	const int pos_full   = name_list_pos("full", nl);
	const bool full_fl   = pos_full  > -1 && nl->inform[pos_full];  // selection with full
	
	const int pos_range  = name_list_pos("range", nl);
	const bool range_fl  = pos_range > -1 && nl->inform[pos_range]; // selection with range
	
	const int pos_slice  = name_list_pos("slice", nl);              // position of slice parameter in select command list
	const bool slice_fl  = pos_slice > -1 && nl->inform[pos_slice]; // selection with slice
	if (slice_fl) slice  = pl->parameters[pos_slice]->double_value; // Parameter has been read. Slice number from select command, if given overwrites slice number which may have been defined by element definition
	
	const int pos_thick  = name_list_pos("thick", nl);              // position of thick flag in select command list
	const bool thick_fl  = pos_slice > -1 && nl->inform[pos_thick]; // selection with slice
	
	if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " i=" << setw(2) << i  << " nl->name=" << nl->name << " full_fl=" << full_fl << " range_fl=" << range_fl
	  << " slice_fl=" << slice_fl << " slice=" << slice << " thick_fl=" << thick_fl << EOL;
	if(full_fl) // use full sequence from start to end, the default
	{
	  nodes[0] = current_sequ->ex_start;
	  nodes[1] = current_sequ->ex_end;
	}
	if(range_fl)
	{
	  if (current_sequ == NULL || current_sequ->ex_start == NULL) // check that there is an active sequence, otherwise crash in get_ex_range
	  {
		warning("makethin range selection without active sequence,", "ignored");
		return 2;
	  }
	  if( get_ex_range(pl->parameters[pos_range]->string, current_sequ, nodes) == 0) // set start nodes[0] and end notes[1] depending on the range string
	  {
		printf("    +++ warning, empty range");
		continue;
	  }
	}
	if(slice_fl || thick_fl) // Set slice number and or thick_fl if present. Add to list of my_selected_elements
	{
	  if(range_fl) // now elements in the sequence in the range
	  {
		node* c_node = nodes[0]; // for range check.  current node
		while (c_node != NULL) // loop over nodes in range,  set slice number in elements
		{
		  element* el_j = c_node->p_elem;
		  const int el_j_slice_pos = name_list_pos("slice",el_j->def->par_names); // position of slice parameter in element list
		  const int el_j_thick_pos = name_list_pos("thick",el_j->def->par_names); // position of thick flag      in element list
		  if (pass_select(el_j->name, slice_select->commands[i]) != 0) // selection on class and pattern done in pass_select. element el_j selected
		  { // the element el_j passes the selection
			if(el_j_slice_pos > -1) el_j->def->par->parameters[el_j_slice_pos]->double_value=slice; // Set the element slice number to the number of slices given in the select statement.
			if(el_j_thick_pos > -1) el_j->def->par->parameters[el_j_thick_pos]->double_value=pl->parameters[pos_thick]->double_value; // Set the element thick flag to what is given in the select statement
		  } // selection
		  if (c_node == nodes[1]) break; // done with last node
		  c_node = c_node->next;
		} // end of while loop over nodes in range
	  } // range_fl
	  else // no range_fl
	  {
		for(int j=0; j< element_list->curr; ++j) // loop over element_list
		{
		  element* el_j = element_list->elem[j];
		  const int el_j_slice_pos = name_list_pos("slice",el_j->def->par_names);
		  const int el_j_thick_pos = name_list_pos("thick",el_j->def->par_names); // position of thick flag      in element list
		  if (pass_select(el_j->name, slice_select->commands[i]) != 0) // selection on class and pattern done in pass_select. element el_j selected
		  { // the element el_j passes the selection
			if(el_j_slice_pos > -1) el_j->def->par->parameters[el_j_slice_pos]->double_value=slice; // Set the element slice number to the number of slices given in the select statement.
			if(el_j_thick_pos > -1) el_j->def->par->parameters[el_j_thick_pos]->double_value=pl->parameters[pos_thick]->double_value; // Set the element thick flag to what is given in the select statement
			if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " el_j->name=" << left << setw(15) << el_j->name << right << " el_j_slice_pos=" << el_j_slice_pos << " el_j_thick_pos=" << el_j_thick_pos << EOL;
		  } // selection
		} // loop over element_list
	  } // range_fl
	} // slice_fl
  } // end of loop over select slice commands
  if (debug_fl()) dump_slices(); // shows where 2 or more slices were selected
  return 0;
}

void makethin(in_cmd* cmd) // public interface to slice sequence, called by exec_command from mad_cmd.c
{
  double CPU_start=clock();
  name_list* nl = cmd->clone->par_names;
  command_parameter_list* pl = cmd->clone->par;
#if __cplusplus >= 201103L // C++11 on
  if (debug_fl()) cout << "using the makethin C++ version with C++11 enabled" << EOL;
#else
  cout << "*** warning *** your compiler/options are not C++11 compatible" << EOL;
#endif
  const int ipos_style = name_list_pos("style", nl);
  string thin_style;
  if (nl->inform[ipos_style] &&  pl->parameters[ipos_style]->string )
  {
	thin_style = buffer(pl->parameters[ipos_style]->string) ;
	cout << "makethin: style chosen : " << thin_style << EOL;
  }

  if(debug_fl() && kill_fringe_fl)   cout << "kill_fringe_fl="   << kill_fringe_fl   << " is on. Flags kill_ent_fringe kill_exi_fringe will be set to true for thick bend body slices" << EOL;
  if(debug_fl() && dipedge_h1_h2_fl) cout << "dipedge_h1_h2_fl=" << dipedge_h1_h2_fl << " is on. Higher order h1, h2 parameters will be kept. Tracking may become non-simplectic" << EOL;

  // first check makethin parameters which influence the selection
  const int ipos_mp = name_list_pos("minimizeparents", nl);
  int MinPar=0;  // or =1  to set minimizeparents to true by default
  if( ipos_mp > -1 && nl->inform[ipos_mp]) MinPar=pl->parameters[ipos_mp]->double_value;
  set_option(const_cast<char*>("minimizeparents"), &MinPar);

  const int ipos_mk = name_list_pos("makeconsistent", nl);
  if( ipos_mk > -1 && nl->inform[ipos_mk])
  {
	int MakeCons=pl->parameters[ipos_mk]->double_value;
	set_option(const_cast<char*>("makeconsistent"), &MakeCons);
  }

  const int ipos_md = name_list_pos("makedipedge", nl);
  if( ipos_md > -1 && nl->inform[ipos_md])
  {
	int iMakeDipedge=pl->parameters[ipos_md]->double_value;
	if (verbose_fl()) cout << "makethin makedipedge flag ipos_md=" << ipos_md << " iMakeDipedge=" << iMakeDipedge << EOL;
	set_option(const_cast<char*>("makedipedge"), &iMakeDipedge);
  }

  if (slice_select->curr > 0)
  {
	int iret=set_selected_elements(); // makethin selection
	if (debug_fl() && iret) cout << "set_selected_elements iret=" << iret << EOL;
  }
  else  warning("makethin: no selection list,","slicing all to one thin lens.");

  if(get_option(const_cast<char*>("makeconsistent"))) force_consistent_slices();

  const int ipos_seq = name_list_pos("sequence", nl);
  char* name = NULL;
  if (nl->inform[ipos_seq] && (name = pl->parameters[ipos_seq]->string))
  {
	const int ipos2 = name_list_pos(name, sequences->list);
	if (ipos2 >= 0)
	{
	  sequence* thick_sequ = sequences->sequs[ipos2];
	  sequence* thin_sequ = slice_sequence(thin_style,thick_sequ); // slice the sequence
	  disable_line(thin_sequ->name, line_list);

	  // 2014-Jul-31  18:13:06  ghislain: make the sequence circular for closed machine by default
	  thin_sequ->start->previous = thin_sequ->end;
	  if (debug_fl()) cout << "makethin: thin_sequ->start->name '" << thin_sequ->start->name << "', thin_sequ->end->name '" << thin_sequ->end->name << "'" << EOL;
	  if (debug_fl()) cout << "makethin: thin_sequ->start->previous->name '" << thin_sequ->start->previous->name << "', thin_sequ->end->next->name '" << thin_sequ->end->next->name << "'" << EOL;
	}
	else warning("unknown sequence ignored:", name);
  }
  else warning("makethin without sequence:", "ignored");
  if (debug_fl()) cout << "makethin: finished in " << (clock()-CPU_start)/CLOCKS_PER_SEC << " seconds" << EOL;
}

//--------  SliceDistPos
SliceDistPos::SliceDistPos(const int n,const bool teapot_fl) : delta(0.5), Delta(0)
{
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
}

void SliceDistPos::Print() const
{
  cout << "SliceDistPos::Print teapot_fl=" << teapot_fl << " n=" << n << " delta=" << delta << " Delta=" << Delta << EOL;
}

//--------  OneElementWithSlices
OneElementWithSlices::OneElementWithSlices(const element* thick_elem,element* thin_elem) // OneElementWithSlices constructor
{
  this->thick_elem=thick_elem;      // for each thick
  sliced_elem.push_back(thin_elem); // there can be several slices
  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " OneElementWithSlices constructor called  thick_elem->name=" <<  thick_elem->name << EOL;
}

//--------  ElementListWithSlices
ElementListWithSlices::ElementListWithSlices(unsigned int verbose) // OneElementWithSlices constructor
{
  // init counters
  get_thin_calls=0;
  get_thin_iteractions=0;
  ilast1=-1;
  ilast2=-1;
  this->verbose=verbose;
  if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " ElementListWithSlices constructor called" << EOL;
}

ElementListWithSlices::~ElementListWithSlices() // destructor
{
  if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " ElementListWithSlices destructor called VecElemWithSlices.size()=" << VecElemWithSlices.size() << EOL;
  for(unsigned int iel=0; iel<VecElemWithSlices.size(); ++iel) delete VecElemWithSlices[iel]; // undo the new OneElementWithSlices(thick_elem,thin_elem);
}

void ElementListWithSlices::PrintCounter() const
{
  cout << "ElementListWithSlices::PrintCounter " << " get_thin_calls=" << get_thin_calls << " get_thin_iteractions=" << get_thin_iteractions;
  if(VecElemWithSlices.size()>0 && get_thin_calls>0) cout << " get_thin_iteractions/get_thin_calls=" << get_thin_iteractions/get_thin_calls << " ineff=" << get_thin_iteractions/((double)VecElemWithSlices.size()*get_thin_calls);
  cout << EOL;
}

int ElementListWithSlices::find_thick(const element* thick_elem) // find thick_element pointer in VecElemWithSlices
{
  get_thin_calls++;
  int ifound=-1;
  if(VecElemWithSlices.size()>0)
  {
	// verbose=3;
	if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " VecElemWithSlices.size()=" << setw(4) << VecElemWithSlices.size()
	  << " " << left << setw(20) << thick_elem->name << right;
	unsigned int isearched=0; // local counter of how many searched
	// look for the element in the list
	if(ilast2 > -1 && VecElemWithSlices[ilast2]->thick_elem==thick_elem)
	{
	  ifound=ilast2; // same as ilast2, no search needed
	  if(verbose>2) cout << " same as ilast2";
	}
	else if(ilast1 > -1 && VecElemWithSlices[ilast1]->thick_elem==thick_elem)
	{
	  ifound=ilast1; // same as ilast1, no search needed
	  if(verbose>2) cout << " same as ilast1";
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
	  if(verbose>2) cout << " not yet known, searched=" << isearched << " get_thin_iteractions=" << get_thin_iteractions << " now ilast1=" << ilast1 << " ilast2=" << ilast2;
	}
	else // found
	{
	  if(verbose>2) cout << " found =" << setw(4) << ifound
		<< " ilast1=" << setw(4) << ilast1 << " ilast2=" << setw(4) << ilast2 << " searched=" << isearched << " get_thin_iteractions=" << get_thin_iteractions;
	  ilast2=ilast1;
	  ilast1=ifound; // remember last found
	}
  }
  return ifound;
}

element* ElementListWithSlices::find_slice(const element* thick_elem,const int slice) // find address of thin slice by slice number for thick_elem
{
  element* result=NULL;

  const int iel=find_thick(thick_elem);
  if(iel<0)
  {
	if(verbose>1) cout << EOL;
	return NULL;
  }
  // thick element found, now check if slice already defined
  int islice=slice-1;
  int nslices=VecElemWithSlices[iel]->sliced_elem.size();

  if(islice < nslices)
  {
	result=VecElemWithSlices[iel]->sliced_elem[islice]; // already done
	if(verbose>1) cout << " result=" << result << " name=" << result->name << EOL;
  }
  else if(verbose>1) cout << " slice " << slice << " still to do" << EOL;
  return result;
}

element* ElementListWithSlices::find_slice(const element* thick_elem,const string& name) // find address of thin slice by slice name for thick_elem
{
  element* result=NULL;

  const int iel=find_thick(thick_elem);
  if(iel<0)
  {
	if(verbose>1) cout << " find_slice=" << left << setw(20) << name << " thick_elem=" << setw(20) << thick_elem->name << right << " not (yet) known" << EOL;
	return NULL;
  }
  const int nslices=VecElemWithSlices[iel]->sliced_elem.size();
  for(unsigned int i=0; i<(unsigned int)nslices; ++i)
  {
	if( string(VecElemWithSlices[iel]->sliced_elem[i]->name)==name) // found
	{
	  result=VecElemWithSlices[iel]->sliced_elem[i]; // can still by NULL, in case of edge elements for e1=0
	  if(verbose>1) cout << " found=" << name << EOL;
	  break;
	}
  }
  if(verbose>1 && result==NULL) cout << " find_slice returns NULL for " << name << EOL;
  return result;
}

void ElementListWithSlices::put_slice(const element* thick_elem,element* thin_elem) // add thin_elem to the list
{
  bool found=false;
  // verbose=3; // CSPE
  if(thick_elem && thin_elem)
  {
    if(verbose>2) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " VecElemWithSlices.size()=" << setw(4) << VecElemWithSlices.size()
      << " thick=" << left << setw(20) << thick_elem->name << setw(20) << " thin=" << thin_elem->name << right << EOL;
    for(unsigned int iel=0; iel<VecElemWithSlices.size(); ++iel)
    {
      if( strcmp(VecElemWithSlices[iel]->thick_elem->name,thick_elem->name) == 0 )
      {
        found=true;
        VecElemWithSlices[iel]->sliced_elem.push_back(thin_elem);
        if(verbose>1) cout << "put_slice found thick name=" << setw(20) << thick_elem->name << " slice name=" << thin_elem->name << " in list at iel=" << iel << " #slices=" << VecElemWithSlices[iel]->sliced_elem.size() << EOL;
        break;
      }
    }
    if(!found)
    {
      OneElementWithSlices* aSliceList= new OneElementWithSlices(thick_elem,thin_elem);
      VecElemWithSlices.push_back(aSliceList);
      if(verbose>1) cout << "put_slice add  thick=" << left << setw(20) << thick_elem->name << setw(20) << " thin=" << thin_elem->name << right << " to list, now VecElemWithSlices.size()=" << VecElemWithSlices.size() << EOL;
    }
  }
  else cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " *** warning *** put_slice called with undefined thick_elem=" << thick_elem << " or thin_elem=" << thin_elem << EOL;
  return;
}

void ElementListWithSlices::Print() const
{
  if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " VecElemWithSlices.size()=" << VecElemWithSlices.size() << EOL;
  cout << "   iel  #slices   base_type         parent_name   parent->base_type    slice_elem->name   slices" << EOL;
  for(unsigned int iel=0; iel<VecElemWithSlices.size(); ++iel) // original
  {
	unsigned int nslices=VecElemWithSlices[iel]->sliced_elem.size();
	if(verbose>1 || nslices>1) // show only if more than 1 slice,  show all in case of verbose
	{
	  const element* el_thick=VecElemWithSlices[iel]->thick_elem;
	  cout << setw(4) << iel << setw(8) << nslices << setw(15) << el_thick->base_type->name << setw(20) << el_thick->name;
	  if(el_thick && el_thick->parent)            cout << setw(20) << el_thick->parent->name; else cout << setw(20) << " ";
	  if(el_thick && el_thick->parent->base_type) cout << setw(20) << el_thick->parent->base_type->name; else cout << setw(20) << " ";	
	  for(unsigned int i=0; i<nslices; ++i)
	  {
		const element* eli=VecElemWithSlices[iel]->sliced_elem[i];
		if(eli) cout << setw(20) << eli->name; else cout << setw(20) << " ";
		cout << " address "  << setw(12) <<  eli;
	  }
	  cout << EOL;
	}
  }
  cout << EOL;
}

//--------  SeqElList
SeqElList::SeqElList(const string& seqname,const string& thin_style,sequence* thick_sequ,sequence* thin_sequ,node* thick_node) : verbose(0),eps(1.e-15),MakeDipedge(true) // SeqElList constructor, eps used to check if values are is compatible with zero
{
  this->seqname=seqname;
  this->thin_style=thin_style;
  this->thick_sequ=thick_sequ;
  this->thin_sequ=thin_sequ;
  this->thick_node=thick_node;
  if(debug_fl())   verbose=1;
  if(verbose_fl()) verbose=2;
  theSliceList=new ElementListWithSlices(verbose);
  theBendEdgeList =new ElementListWithSlices(verbose);
  // verbose=3; // -- special for code development ---   turn on extra debugging within SeqElList
  MakeDipedge=get_option(const_cast<char*>("makedipedge"));
  if(verbose && !MakeDipedge) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " ***  makedipedge is off, should always be on except for backwards compatibility checks  *** " << EOL;
  if(verbose>1)
  {
	cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " constructor seqname=" << seqname << " makedipedge check  MakeDipedge=" << MakeDipedge << EOL;
  }
}

SeqElList::~SeqElList() // destructor
{
  delete theSliceList;
  delete theBendEdgeList;
}

double SeqElList::simple_at_shift(const int slices,const int slice_no) // return at relative shifts from centre of unsliced magnet
{
  const int n = slices;
  const int i = slice_no;
  return n>1 ? (2.0*i-1)/(2.0*n)-0.5 : 0.0;
}

double SeqElList::teapot_at_shift(int slices, int slice_no)
{
  const int n = slices;
  const int i = slice_no;
  return n>1 ? 0.5*n*(1-2*i+n)/(1.0-n*n) : 0.0; // see  http://ab-dep-abp.web.cern.ch/ab-dep-abp/LCU/LCU_meetings/2012/20120918/LCU_makethin_2012_09_18.pdf
}

double SeqElList::collim_at_shift(const int slices,const int slice_no)
{
  const int n = slices;
  const int i = slice_no;
  return n>1 ? (i-1.0)/(n-1.0)-0.5 : 0.0;
}

double SeqElList::default_at_shift(const int slices, const int slice_no)
{
  return slices>4 ? simple_at_shift(slices, slice_no) : teapot_at_shift(slices, slice_no);
}

double SeqElList::at_shift(const int slices,const int slice_no,const string& local_thin_style) // return relative shifts from centre of unsliced magnet
{
  if(verbose_fl() && local_thin_style!=thin_style) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " local_thin_style=" << local_thin_style << " thin_style=" << thin_style << EOL;

  if (!slices || !slice_no) fatal_error("makethin: invalid slicing for zero slices",local_thin_style.c_str());
  if      (local_thin_style==string(""))        return default_at_shift(slices,slice_no);
  else if (local_thin_style==string("simple"))  return simple_at_shift(slices,slice_no);
  else if (local_thin_style==string("teapot"))  return teapot_at_shift(slices,slice_no);
  else if (local_thin_style==string("collim"))  return collim_at_shift(slices,slice_no);
  else fatal_error("makethin: Style chosen not known:",local_thin_style.c_str());
  return 0;
}

void SeqElList::Print() const
{
  cout << "SeqElList::Print seqname=" << seqname << " theSliceList->VecElemWithSlices.size()=" << theSliceList->VecElemWithSlices.size() << " thin_style=\"" << thin_style << "\"" << EOL;

  cout << EOL << "   theSliceList:" << EOL; theSliceList->Print();
  if(verbose) theSliceList->PrintCounter();

  cout << EOL << "theBendEdgeList:" << EOL; theBendEdgeList->Print();
  if(verbose) theBendEdgeList->PrintCounter();
}

int SeqElList::translate_k(command_parameter* *kparam, command_parameter* *ksparam,const command_parameter* angle_param, command_parameter* kn_param, command_parameter* ks_param)
// translate k0,k1,k2,k3 & k0s,k1s,k2s,k3s to kn{} and ks{}
{
  int angle_conversion=0;

  if ((kparam == NULL) && (ksparam == NULL)) fatal_error("translate_k: no kparams to convert","");

  if (angle_param) // if we have an angle we ignore any given k0
  {
	kparam[0] =  new_command_parameter(const_cast<char*>("k0"), k_double);
	angle_conversion=1; // note we do not divide by length, just to multiply again afterwards
	if (angle_param->expr) kparam[0]->expr =  clone_expression(angle_param->expr);
	kparam[0]->double_value = angle_param->double_value;
  }

  for(int i=0; i<4; ++i) // initialize the parameters with NULL for expression and 0 for the value
  {
	kn_param->expr_list->list[i] = NULL; kn_param->double_array->a[i] = 0;
	ks_param->expr_list->list[i] = NULL; ks_param->double_array->a[i] = 0;
	if (kparam[i]) // copy existing k's
	{
	  if (kparam[i]->expr) kn_param->expr_list->list[i] = clone_expression(kparam[i]->expr);
	  kn_param->double_array->a[i] = kparam[i]->double_value;
	}
	if (ksparam[i])
	{
	  if (ksparam[i]->expr) ks_param->expr_list->list[i] = clone_expression(ksparam[i]->expr);
	  ks_param->double_array->a[i] = ksparam[i]->double_value;
	}
	kn_param->expr_list->curr++; kn_param->double_array->curr++; // update the number of k's in our arrays
	ks_param->expr_list->curr++; ks_param->double_array->curr++;
  }
  return angle_conversion;
}

element* SeqElList::create_thick_slice(element* thick_elem,const int i) // create entry/exit slice _s,  and body _b elements like   mqxa.1r1_s: quadrupole, l:=l.mqxa   * 1/10, k1:=kqx.r1 + ktqx1.r1;      and add to the sequence
{
  const int n = get_slices_from_elem(thick_elem);
  const char* eltype=thick_elem->base_type->name;
  string slice_name;
  const bool entry_fl=(i==0);
  const bool  exit_fl=(i==n);
  const bool IsBend=(string(eltype).substr(1,4)==string("bend"));
  const bool IsRbend= string(thick_elem->base_type->name) == string("rbend");
  if      (entry_fl)  slice_name=string(thick_elem->name)+"_en"; // entry
  else if (exit_fl)   slice_name=string(thick_elem->name)+"_ex"; // exit
  else                slice_name=string(thick_elem->name)+"_bo"; // body  slice

  element* slice_elem;
  if( (slice_elem = theSliceList->find_slice(thick_elem,slice_name)))
  {
    if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " slice_name already exists, use it" << EOL;
    return slice_elem;
  }

  if(IsRbend)
  {
    cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << eltype << " create " << slice_name << " based on " << thick_elem->name << " *** warning *** slicing of thick rbends requested -- generates extra fringe fields from the angle at the surface " << EOL;
  }

  SliceDistPos SP(n, thin_style==string("teapot") );
  if(verbose>1)
  {
    cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << eltype << " create " << slice_name << " based on " << thick_elem->name
    << " i=" << i << " n=" << n
    << " entry_fl=" << entry_fl
    << " exit_fl="  << exit_fl
    << " IsBend=" << IsBend
    << " IsRbend=" << IsRbend
    << " " << my_dump_element(thick_elem); SP.Print();
    if(thick_elem->parent) cout << " dump also parent " << thick_elem->parent->name << my_dump_element(thick_elem->parent);
  }

  expression* l_par_expr=my_get_param_expression(thick_elem,const_cast<char*>("l")); // get length expression
  expression* angle_par_expr = my_get_param_expression(thick_elem,const_cast<char*>("angle")); // get angle expressions - relevant for bends

  if(l_par_expr==NULL) // compound_expr with scaling will fail   -- should never happen
  {
    cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " *** error *** l_par_expr=" << l_par_expr << EOL;
    exit(1);
  }

  double LengthFraction;
  if(entry_fl || exit_fl) LengthFraction=SP.delta; // start / end slice
  else                    LengthFraction=SP.Delta; // the middle or body piece

  l_par_expr                        = compound_expr(l_par_expr,    0,const_cast<char*>("*"), NULL,LengthFraction); // multiply length parameter expression with LengthFraction
  if(angle_par_expr) angle_par_expr = compound_expr(angle_par_expr,0,const_cast<char*>("*"), NULL,LengthFraction); // multiply angle  parameter expression with LengthFraction, only relevant for bends

  command* cmd = clone_command(thick_elem->def); 	// clone existing command to define the new element, result something like  mqxa.1r1_s: quadrupole,polarity:= 1,k1:=kqx.r1 + ktqx1.r1 ;

  if(verbose>1)
  {
    cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " create_thick_slice after clone for " << thick_elem->name << " ";
    SP.Print();
    cout  << " LengthFraction=" << LengthFraction << " scaled l_par_expr " << my_dump_expression(l_par_expr);
    if(angle_par_expr) cout << " scaled angle_par_expr " << my_dump_expression(angle_par_expr);
    cout << EOL;
  }

  if( IsRbend && rbarc_fl() ) // straight length rbend, needs translation
  {
    double angle = expression_value(angle_par_expr, 2);
    double factor = 2*sin(angle/2) / angle; // <= 1
    double val=LengthFraction*thick_elem->length*factor; // new reduced straight rbend slice length
    char** toks = tmp_l_array->p;
    ostringstream ostr;
    ostr << setprecision(15) << val; // use the value as string
    toks[0]=CopToNewC_String(ostr.str());
    l_par_expr = make_expression(1,toks); // new simple expression
    if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " thick rbend slicing, scale thick curved length " << thick_elem->length << " by factor=" << setprecision(15) << factor << " to val=" << val << EOL;
  }

  // now use the scaled length and if relevant angle parameters to set up the new slice_elem via cmd
  const int length_i = name_list_pos(const_cast<char*>("l"),cmd->par_names);
  if(length_i > -1)
  {
    cmd->par->parameters[length_i]->expr=l_par_expr; // use the length expression in cmd to set up slice_elem
  }
  else
  {
    cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " *** error *** thick_elem " <<  thick_elem->name << " has no length parameter : " << my_dump_element(thick_elem);
    return NULL;
  }
  const int angle_i = name_list_pos(const_cast<char*>("angle"),cmd->par_names);
  if(verbose>1) cout << EOL << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " angle_i=" << angle_i << EOL;
  if(angle_i > -1)
  {
    cmd->par->parameters[angle_i]->expr=angle_par_expr; // use the scaled angle_par_expr in cmd to set up slice_elem
  }

  if (verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " done thin_style=\"" << thin_style << "\"" << " i=" << i << " new cmd for slice_elem : " << my_dump_command(cmd);

  slice_elem = make_element(const_cast<char*>(slice_name.c_str()),const_cast<char*>(eltype),cmd,-1); // make new element (without parent) using the command cmd, -1 means avoid warnings
  if (verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ <<  my_dump_element(slice_elem) << EOL;

  ParametersActiveOn(slice_elem); // Activate attributes, important when derived from parents -- turns also slice number on - not wanted

  ParameterRemove(const_cast<char*>("slice"),slice_elem); // slicing done, no reason to leave the slice parameter

  SetParameterValue("thick",slice_elem,true,k_logical);
  ParameterTurnOn("thick",slice_elem); //-- so that thick=true is written  to signal this for thick tracking
  if (verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " done create_thick_slice slice_name=" << slice_name << " from thick element " << thick_elem->name << " n=" << n  << " : " << my_dump_element(slice_elem) << EOL;

  if(IsBend)
  {
    if(entry_fl) // bend entry,  remove/turn off exit parameters
    {
      ParameterRemove("e2"    ,slice_elem);
      ParameterRemove("h2"    ,slice_elem);
      ParameterTurnOn("fint" ,slice_elem);
      ParameterTurnOn("fintx",slice_elem);
      SetParameterValue("fintx",slice_elem,0); // leave with value 0, otherwise taking fint
    }
    else if(exit_fl) // bend exit, remove entry parameters
    {
      ParameterRemove(const_cast<char*>("e1")   ,slice_elem);
      int i_fint = name_list_pos(const_cast<char*>("fint"),slice_elem->def->par_names);
      const bool fint_on =  (i_fint > -1) && slice_elem->def->par_names->inform[i_fint];
      int i_fintx = name_list_pos(const_cast<char*>("fintx"),slice_elem->def->par_names);
      const bool fintx_on =  (i_fintx > -1) && slice_elem->def->par_names->inform[i_fintx];
      if(fintx_on) ParameterRemove("fint",slice_elem); // fintx is on and will be used, just remove any fint on the exit
      else if(fint_on) // not fintx, use fint as fintx for exit
      {
        ParameterTurnOn("fintx",slice_elem);
        if(i_fintx) // should be there, just inform off
        {
          bool found=false;
          double fint_value=my_get_int_or_double_value(slice_elem,"fint",found);
          SetParameterValue("fintx",slice_elem,fint_value);
          if (verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " no fintx, use fint value " << fint_value << " as fintx for exit" << EOL;
        }
        ParameterRemove("fint",slice_elem); // remove fint on exit
      }
      ParameterRemove("h1",slice_elem);
    }
    else Remove_All_Fringe_Field_Parameters(slice_elem); // thick magnet body, remove fringe fields
  }
  if (verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ <<  my_dump_element(slice_elem) << EOL;
  theSliceList->put_slice(thick_elem,slice_elem); //-- store what is done in theSliceList
  return slice_elem;
} // create_thick_slice

element* SeqElList::create_sliced_magnet(const element* thick_elem, int slice_no,bool ThickSLice) // create the sliced element, recursively for parents
{
  element *thin_elem_parent;
  int slices = get_slices_from_elem(thick_elem);
  if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << thick_elem->name << " slices=" << slices << " ThickSLice=" << ThickSLice << " slice_no=" << slice_no << EOL;

  if (thick_elem == thick_elem->parent) return NULL; // no further parent to consider
  else
  {
    if(verbose>1) printf("recursively slice parent:");
    thin_elem_parent = create_sliced_magnet(thick_elem->parent,slice_no,ThickSLice); // recursively slice parent
  }
  const command_parameter* at_param = return_param(("at"),thick_elem); 	// handle parent with possibly different slice number than child
  const int minimizefl=get_option(const_cast<char*>("minimizeparents")) && !at_param && thick_elem == thick_elem->parent;
  if(minimizefl)
  {
    slice_no=slices=1; // do not slice this one
  }
  if(slice_no > slices && thick_elem!=thick_elem->parent ) slice_no=1; // check, but not for base classes
  if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " creating new kn_param, ks_param expr_list" << EOL;

  element *thin_elem;
  if( (thin_elem = theSliceList->find_slice(thick_elem,slice_no)) ) return thin_elem; // already done
  const command_parameter* length_param = return_param_recurse(("l"),thick_elem);
  const command_parameter* angle_param  = return_param_recurse(("angle"),thick_elem);
  command_parameter* kparam[4]; // clone to allow for changes in copy, use my clone version which does not crash if called with NULL
  kparam[0]    = clone_command_parameter(return_param_recurse(("k0"),thick_elem) );
  kparam[1]    = clone_command_parameter(return_param_recurse(("k1"),thick_elem) );
  kparam[2]    = clone_command_parameter(return_param_recurse(("k2"),thick_elem) );
  kparam[3]    = clone_command_parameter(return_param_recurse(("k3"),thick_elem) );
  command_parameter* ksparam[4];
  ksparam[0]   = clone_command_parameter(return_param_recurse(("k0s"),thick_elem) );
  ksparam[1]   = clone_command_parameter(return_param_recurse(("k1s"),thick_elem) );
  ksparam[2]   = clone_command_parameter(return_param_recurse(("k2s"),thick_elem) );
  ksparam[3]   = clone_command_parameter(return_param_recurse(("k3s"),thick_elem) );
  command_parameter* kn_param = clone_command_parameter( return_param_recurse(("knl"),thick_elem) );
  command_parameter* ks_param = clone_command_parameter( return_param_recurse(("ksl"),thick_elem) );

  int knl_flag = 0,ksl_flag = 0;
  if (kn_param) {kn_param = clone_command_parameter(kn_param); knl_flag++;}
  if (ks_param) {ks_param = clone_command_parameter(ks_param); ksl_flag++;}

  int angle_conversion = 0;

  if ((kparam[0] || kparam[1] || kparam[2] || kparam[3] || angle_param || ksparam[0] || ksparam[1] || ksparam[2] || ksparam[3])
      && (kn_param==NULL && ks_param==NULL)) 	// translate k0,k1,k2,k3,angles
  {
    if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " creating new kn_param, ks_param expr_list" << EOL;
    kn_param = new_command_parameter(const_cast<char*>("knl"), k_double_array);
    kn_param->expr_list = new_expr_list(10);
    kn_param->double_array = new_double_array(10);
    ks_param = new_command_parameter(const_cast<char*>("ksl"), k_double_array);
    ks_param->expr_list = new_expr_list(10);
    ks_param->double_array = new_double_array(10);
    angle_conversion = translate_k(kparam,ksparam,angle_param,kn_param,ks_param);
  }

  kn_param = scale_and_slice(kn_param,length_param,slices,angle_conversion,knl_flag+ksl_flag);
  ks_param = scale_and_slice(ks_param,length_param,slices,angle_conversion,knl_flag+ksl_flag);
  if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << thick_elem->name << EOL;

  if(ThickSLice)
  { // from  knl:={amb ,( kmb ) * ( lmb ) };   to
    if(verbose>1) printf("debug %s %s line %d ThickSLice=%d set kn, ks to zero\n",__FILE__,__FUNCTION__,__LINE__,ThickSLice);
    kn_param=NULL;
    ks_param=NULL;
  }

  // set up new multipole command
  command* cmd = new_command( CopToNewC_String("thin_multipole"), 20, 20, CopToNewC_String("element"), CopToNewC_String("none"), 0, 8); // 0 is link, multipole is 8, see "multipole: element none 0 8 " in mad_dict.c
  add_cmd_parameter_new(cmd,1.,"magnet",0); // parameter magnet with value of 1 and inf=0         not really needed ?   is the default, see mad_dict.c
  if(!minimizefl)
  {
    add_cmd_parameter_clone(cmd,return_param("at"  ,thick_elem),"at"  ,1);
    add_cmd_parameter_clone(cmd,return_param("from",thick_elem),"from",1);
    add_lrad(cmd,length_param,slices); // add l, lrad
    add_cmd_parameter_clone(cmd,(const command_parameter*)kn_param,"knl",1);
    add_cmd_parameter_clone(cmd,(const command_parameter*)ks_param,"ksl",1);
  }
  if(verbose>1) cout << my_dump_command(cmd) << EOL; // magnet, l, lrad defined
  //--- now the arguments which are copied from the thick element
  add_cmd_parameter_clone(cmd,return_param_recurse(("apertype"),thick_elem),("apertype"),1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("aperture"),thick_elem),("aperture"),1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("aper_tol"),thick_elem),("aper_tol"),1);
  add_cmd_parameter_clone(cmd,return_param(        ("bv"),      thick_elem),("bv"),      1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("tilt"),    thick_elem),("tilt"),    1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("kmax"),    thick_elem),("kmax"),    1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("kmin"),    thick_elem),("kmin"),    1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("calib"),   thick_elem),("calib"),   1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("polarity"),thick_elem),("polarity"),1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("mech_sep"),thick_elem),("mech_sep"),1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("v_pos"),   thick_elem),("v_pos"),   1);
  if(verbose>1) cout << my_dump_command(cmd) << EOL;
  // create element with this command
  char* thin_name;
  if (slices==1 && slice_no==1) thin_name=buffer(CopToNewC_String(thick_elem->name));
  else
  {
    thin_name = make_thin_name(thick_elem->name,slice_no);
    if(verbose>1) printf("verbose %s %s line %d make_thin_name(%s,%d)=%s\n",__FILE__,__FUNCTION__,__LINE__,thick_elem->name,slice_no,thin_name);
  }
  if (thin_elem_parent)
  {
    if(verbose>1) printf("verbose %s %s line %d make_element(%s,%s,cmd,-1);\n",__FILE__,__FUNCTION__,__LINE__,thin_name,thin_elem_parent->name);
    thin_elem = make_element(thin_name,thin_elem_parent->name,cmd,-1);
  }
  else
  {
    if(verbose>1) printf("verbose %s %s line %d make_element(%s,\"multipole\",cmd,-1);\n",__FILE__,__FUNCTION__,__LINE__,thin_name);
    thin_elem = make_element(thin_name,const_cast<char*>("multipole"),cmd,-1);
  }
  thin_elem->length = 0;
  thin_elem->bv = el_par_value(const_cast<char*>("bv"),thin_elem);
  if (thin_elem_parent && thin_elem_parent->bv)
  {
    thin_elem->bv = thin_elem_parent->bv;
  }
  if(verbose>1) printf("verbose %s %s line %d put_slice(thick_elem %s, thin_elem %s,%d);\n",__FILE__,__FUNCTION__,__LINE__,thick_elem->name,thin_elem->name,slice_no);
  theSliceList->put_slice(thick_elem,thin_elem);
  return thin_elem;
}

element* SeqElList::create_thin_solenoid(const element* thick_elem, int slice_no) // create thin solenoid element, similar to create_sliced_magnet
{
  element *thin_elem_parent;
  if (thick_elem == thick_elem->parent) return NULL;
  else thin_elem_parent = create_thin_solenoid(thick_elem->parent,slice_no); // recursively slice parent
  element *thin_elem;
  if((thin_elem = theSliceList->find_slice(thick_elem,slice_no))) return thin_elem; // check to see if we've already done this one
  // get parameters from the thick solenoid element
  const command_parameter* length_param  = return_param_recurse(("l"),thick_elem);
  const command_parameter* ks_param      = return_param_recurse(("ks"),thick_elem);
  const command_parameter* at_param      = return_param(("at"),thick_elem);

  int slices = get_slices_from_elem(thick_elem);
  const int minimizefl=get_option(const_cast<char*>("minimizeparents")) && !at_param && thick_elem == thick_elem->parent;
  if(minimizefl) slice_no=slices=1; // do not slice this one

  // set up new solenoid command
  command* cmd = new_command(buffer(const_cast<char*>("thin_solenoid")), 20, 20, // max num names, max num param
                             buffer(const_cast<char*>("element")), buffer(const_cast<char*>("none")), 0, 9); // 0 is link, solenoid is 9
  add_cmd_parameter_new(cmd,1.,"magnet",0); // parameter magnet with value of 1 and inf=0

  if(!minimizefl)
  {
	add_cmd_parameter_clone(cmd,return_param(("at")  ,thick_elem),("at")  ,1);
	add_cmd_parameter_clone(cmd,return_param(("from"),thick_elem),("from"),1);
	add_lrad(cmd,length_param,slices);
  }
  add_cmd_parameter_clone(cmd,ks_param,("ks"),1); // keep ks
  if(!minimizefl)
  {
	if (length_param && ks_param) /* in addition provide   ksi = ks * l /slices */
	{
	  command_parameter* ks_par = cmd->par->parameters[cmd->par->curr] = clone_command_parameter(ks_param); // start from clone of ks
	  strcpy(ks_par->name,"ksi"); // change name to ksi
	  if (length_param->expr || ks_par->expr) // first step is ks * l calculation, expression or value
	  {
		ks_par->expr = compound_expr(ks_par->expr,ks_par->double_value,const_cast<char*>("*"),length_param->expr,length_param->double_value); // multiply expression with length
	  }
	  else ks_par->double_value *= length_param->double_value; // multiply value with length
	  if (slices > 1) // 2nd step, divide by slices, expression or number
	  {
		if (ks_par->expr) ks_par->expr = compound_expr(ks_par->expr,0.,const_cast<char*>("/"),NULL,slices);
		else ks_par->double_value /= slices;
	  }
	  add_to_name_list(const_cast<char*>("ksi"),1,cmd->par_names);
	  cmd->par->curr++;
	}
  }
  add_cmd_parameter_clone(cmd,return_param_recurse(("apertype"),thick_elem),("apertype"),1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("aperture"),thick_elem),("aperture"),1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("aper_tol"),thick_elem),("aper_tol"),1);
  add_cmd_parameter_clone(cmd,return_param(("bv"),thick_elem),("bv"),1);
  add_cmd_parameter_clone(cmd,return_param(("tilt"),thick_elem),("tilt"),1);
  // create element with this command
  char* thin_name;
  if (slices==1 && slice_no==1) thin_name=buffer(CopToNewC_String(thick_elem->name)); //old buffer(thick_elem->name);
  else thin_name = make_thin_name(thick_elem->name,slice_no);
  if (thin_elem_parent)
  {
	thin_elem = make_element(thin_name,thin_elem_parent->name,cmd,-1);
  }
  else
  {
	thin_elem = make_element(thin_name,const_cast<char*>("solenoid"),cmd,-1);
  }
  thin_elem->length = 0;
  thin_elem->bv = el_par_value(const_cast<char*>("bv"),thin_elem);
  if (thin_elem_parent && thin_elem_parent->bv)
  {
	thin_elem->bv = thin_elem_parent->bv;
  }
  theSliceList->put_slice(thick_elem,thin_elem);
  return thin_elem;
}

element* SeqElList::create_thin_elseparator(const element* thick_elem, int slice_no) // create thin elseparator element, similar to create_thin_solenoid
{
  element *thin_elem_parent;
  if (thick_elem == thick_elem->parent) return NULL;
  else thin_elem_parent = create_thin_elseparator(thick_elem->parent,slice_no); // recursively slice parent
  element *thin_elem;
  if((thin_elem = theSliceList->find_slice(thick_elem,slice_no))) return thin_elem; // check to see if we've already done this one

  // get parameters from the thick elseparator element
  const command_parameter* length_param  = return_param_recurse(("l"),thick_elem);
  const command_parameter* ex_param      = return_param_recurse(("ex"),thick_elem);
  const command_parameter* ey_param      = return_param_recurse(("ey"),thick_elem);
  const command_parameter* tilt_param    = return_param_recurse(("tilt"),thick_elem);
  const command_parameter* at_param      = return_param(("at"),thick_elem);

  int slices = get_slices_from_elem(thick_elem);
  const int minimizefl=get_option(const_cast<char*>("minimizeparents")) && !at_param && thick_elem == thick_elem->parent;
  if(minimizefl)
  {
	slice_no=slices=1; /* do not slice this one */
  }

  // set up new elseparator command
  command* cmd = new_command(buffer(const_cast<char*>("thin_elseparator")), 20, 20, // max num names, max num param
                             buffer(const_cast<char*>("element")), buffer(const_cast<char*>("none")), 0, 11); // 0 is link, elseparator is 11
  add_cmd_parameter_new(cmd,1.,"magnet",0); // parameter magnet with value of 1 and inf=0

  if(!minimizefl)
  {
	add_cmd_parameter_clone(cmd,return_param(("at")  ,thick_elem),("at")  ,1);
	add_cmd_parameter_clone(cmd,return_param(("from"),thick_elem),("from"),1);
	add_lrad(cmd,length_param,slices);
  }
  add_cmd_parameter_clone(cmd,ex_param,("ex"),1); // keep ex
  add_cmd_parameter_clone(cmd,ey_param,("ey"),1); // keep ey
  add_cmd_parameter_clone(cmd,tilt_param,("tilt"),1); // keep tilt
  if(!minimizefl)
  { // create ex_l from ex
	if (length_param && ex_param) // in addition provide ex_l = ex * l /slices
	{
	  command_parameter* ex_par = cmd->par->parameters[cmd->par->curr] = clone_command_parameter(ex_param); // start from clone of ex
	  strcpy(ex_par->name,"ex_l"); // change name to ex_l
	  if (length_param->expr && ex_par->expr) // first step is ex * l calculation, expression or value
	  {
		ex_par->expr = compound_expr(ex_par->expr,ex_par->double_value,const_cast<char*>("*"),length_param->expr,length_param->double_value); /* multiply expression with length */
	  }
	  else ex_par->double_value *= length_param->double_value; // multiply value with length
	  if (slices > 1) // 2nd step, divide by slices, expression or number
	  {
		if (ex_par->expr) ex_par->expr = compound_expr(ex_par->expr,0.,const_cast<char*>("/"),NULL,slices);
		else ex_par->double_value /= slices;
	  }
	  add_to_name_list(const_cast<char*>("ex_l"),1,cmd->par_names);
	  cmd->par->curr++;
	}
	// create ey_l from ey
	if (length_param && ey_param) // in addition provide   ey_l = ey * l /slices
	{
	  command_parameter* ey_par = cmd->par->parameters[cmd->par->curr] = clone_command_parameter(ey_param); // start from clone of ey
	  strcpy(ey_par->name,"ey_l"); // change name to ey_l
	  if (length_param->expr && ey_par->expr) // first step is ey * l calculation, expression or value
	  {
		ey_par->expr = compound_expr(ey_par->expr,ey_par->double_value,const_cast<char*>("*"),length_param->expr,length_param->double_value); // multiply expression with length
	  }
	  else ey_par->double_value *= length_param->double_value; // multiply value with length
	  if (slices > 1) // 2nd step, divide by slices, expression or number
	  {
		if (ey_par->expr) ey_par->expr = compound_expr(ey_par->expr,0.,const_cast<char*>("/"),NULL,slices);
		else ey_par->double_value /= slices;
	  }
	  add_to_name_list(const_cast<char*>("ey_l"),1,cmd->par_names);
	  cmd->par->curr++;
	}
  }
  add_cmd_parameter_clone(cmd,return_param_recurse(("apertype"),thick_elem),("apertype"),1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("aperture"),thick_elem),("aperture"),1);
  add_cmd_parameter_clone(cmd,return_param_recurse(("aper_tol"),thick_elem),("aper_tol"),1);
  add_cmd_parameter_clone(cmd,return_param(("bv"),thick_elem),("bv"),1);
  add_cmd_parameter_clone(cmd,return_param(("tilt"),thick_elem),("tilt"),1);
  // create element with this command
  char* thin_name;
  if (slices==1 && slice_no==1) thin_name=buffer(CopToNewC_String(thick_elem->name));
  else thin_name = make_thin_name(thick_elem->name,slice_no);
  if (thin_elem_parent) thin_elem = make_element(thin_name,thin_elem_parent->name,cmd,-1);
  else thin_elem = make_element(thin_name,const_cast<char*>("elseparator"),cmd,-1);
  thin_elem->length = 0;
  thin_elem->bv = el_par_value(const_cast<char*>("bv"),thin_elem);
  if (thin_elem_parent && thin_elem_parent->bv) thin_elem->bv = thin_elem_parent->bv;
  theSliceList->put_slice(thick_elem,thin_elem);
  return thin_elem;
}

void SeqElList::slice_this_node() // main stearing what to do.   called in loop over nodes which can be sliced, makes slices, adds them to  thin_sequ
{
  bool UseDipedges=true; // Normally true.   For tests set false, to see the result without dipedges, dipedges will then be created but not written

  element* thick_elem=thick_node->p_elem; // work directly on this element  --  to do translaton only once, element maybe used several times in sequence
  // element* thick_elem=clone_element(thick_node->p_elem); // work on clone of this element, --  no

  int nslices=get_slices_from_elem(thick_elem);
  if(nslices<1)
  {
	if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << thick_elem->name << " nslices=" << nslices << " is less than one, place thick_node " << thick_node->name << " without any change" << EOL;
	add_node_at_end_of_sequence(thick_node,thin_sequ); // straight copy
	return;
  }

  // verbose=3; // CSPE enable to get extra debug info
  const bool IsQuad =   string(thick_node->base_name) == string("quadrupole");
  const bool IsRbend=   string(thick_node->base_name) == string("rbend");
  const bool IsBend = ( string(thick_node->base_name) == string("sbend") || IsRbend );
  const bool ThickSLice=thick_fl(thick_elem);

  if(verbose>1)
  {
	cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << thick_elem->name << " " << thick_node->base_name << " thin_style=\"" << thin_style << "\""
	<< " nslices=" << nslices << " IsQuad=" << IsQuad << " IsBend=" << IsBend << " ThickSLice=" << ThickSLice << " at_value=" << setw(10) << thick_node->at_value << " from_name=";
	if(thick_node->from_name) cout << thick_node->from_name; else cout << "NULL ";
	cout << EOL << my_dump_element(thick_elem) << EOL;
  }

  if(ThickSLice && nslices>1 && !IsQuad && !IsBend)
  {
	if(nslices>1) cout << "++++++ warning: " << thick_elem->name << " is a " << thick_node->base_name << " nslices=" << nslices << " thick slicing with nslices>1 for " << thick_node->base_name << ". Set nslices=1" << EOL;
    nslices=1;
  }

  element *EntryDipedge=NULL, *ExitDipedge=NULL, *sbend_el=NULL;
  element *en = NULL, *bo = NULL, *ex = NULL; // pointers to thick  entry, body, exit

  if( IsBend && MakeDipedge) // find any existing EntryDipedge, sbend_el, ExitDipedge    and use them
  {	
	EntryDipedge=theBendEdgeList->find_slice(thick_elem,string(thick_elem->name)+"_den"); // dipedge entry, NULL if not yet known or e1=0
	ExitDipedge =theBendEdgeList->find_slice(thick_elem,string(thick_elem->name)+"_dex"); // dipedge exit,  BULL if not yet known or e2=0
	if(IsRbend) sbend_el=theBendEdgeList->find_slice(thick_elem,string(thick_elem->name)+"_sbend");
    if(verbose>1)
    {
      if(verbose>1)    cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__  << "              " << setw(20) <<   thick_elem->name << " " << thick_node->base_name << EOL;
      if(EntryDipedge) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__  << " EntryDipedge=" << setw(20) << EntryDipedge->name << " already exists " << EntryDipedge << EOL;
      if(ExitDipedge)  cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__  << "  ExitDipedge=" << setw(20) <<  ExitDipedge->name << " already exists " << EntryDipedge << EOL;
      if(sbend_el)     cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__  << "     sbend_el=" << setw(20) <<     sbend_el->name << " already exists " << EntryDipedge << EOL;
    }
  }

  if(IsRbend && ThickSLice && sbend_el==NULL && (rbend_to_sbend_fl || MakeDipedge) ) // create new thick sbend_el for this rbend,  translate all thick rbends
  {
    if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__  << " sbend_el=" << string(thick_elem->name)+"_sbend" << " does not (yet) exist, make it"<< EOL;
    sbend_el=sbend_from_rbend(thick_elem,MakeDipedge); // sbend_el for rbend not yet existing, make the element and put it in global element_list
    theBendEdgeList->put_slice(thick_elem,sbend_el); // to remember this has been translated
    if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__  << " sbend_el created sbend_el=" << sbend_el << EOL;
    if(verbose>2) dump_slices();
  } // ThickSLice && IsRbend

  if( IsBend && MakeDipedge && EntryDipedge==NULL) // create new EntryDipedge for this bend
  { // first look if e1 or h1 are there
    const command_parameter   *e1param = return_param(("e1")  ,(const element*) thick_elem);
    const command_parameter   *h1param = return_param(("h1")  ,(const element*) thick_elem);
    const command_parameter *fintparam = return_param(("fint"),(const element*) thick_elem);
    if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " e1param=" << e1param << " cmd_par_val(e1param)=" << cmd_par_val(e1param) << EOL;
    if((fabs(cmd_par_val(e1param))>eps) || (fabs(cmd_par_val(h1param))>eps) || (fabs(cmd_par_val(fintparam))>eps) || IsRbend) // has entrance fringe fields
    {
      EntryDipedge=create_bend_dipedge_element(thick_elem,true); // make new StartEdge element and remove e1 from thick_elem, change rbend to sbend
      if(IsRbend) add_half_angle_to(thick_elem,EntryDipedge,"e1");
      theBendEdgeList->put_slice(thick_elem,EntryDipedge);       // to remember this has been translated
      if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__  << " now EntryDipedge=" << EntryDipedge << EOL;
	}
  }

  if( IsBend && MakeDipedge && ExitDipedge==NULL) // create new ExitDipedge for this bend
  { // first look if e2 or h2 are there
    const command_parameter *e2param    = return_param(("e2"),(const element*) thick_elem);
    const command_parameter *h2param    = return_param(("h2"),(const element*) thick_elem);
    const command_parameter *fintparam  = return_param(("fint"),(const element*) thick_elem);
    const command_parameter *fintxparam = return_param(("fintx"),(const element*) thick_elem);
    if((fabs(cmd_par_val(e2param))>eps) || (fabs(cmd_par_val(h2param))>eps) || (fabs(cmd_par_val(fintparam))>eps) || (cmd_par_val(fintxparam)>eps) || IsRbend)  // has exit fringe fields
    {
      ExitDipedge=create_bend_dipedge_element(thick_elem,false); // make new ExitDipedge element and remove e2 from thick_elem,  change rbend to sbend
      if(IsRbend) add_half_angle_to(thick_elem,ExitDipedge,"e1");
      theBendEdgeList->put_slice(thick_elem,ExitDipedge);   // to remember this has been translated
      if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__  << " now  ExitDipedge=" << EntryDipedge << " " << my_dump_element(ExitDipedge) << EOL;
    }
  } // new ExitDipedge

  if(sbend_el) thick_elem=sbend_el; // in case of rbend -> sbend translation, continue with the translated sbend
  if(EntryDipedge || ExitDipedge) Remove_All_Fringe_Field_Parameters(thick_elem); // remove from body what is now taken care of by dipedges

  // prepare for slicing
  element *sliced_elem=NULL;                  // pointer to new sliced element

  string local_thin_style=thin_style; // work here with a copy of thin_style that can be modified for collim
  if (strstr(thick_node->base_name,"collimator"))
  {
    local_thin_style = "collim"; // for collimators change slice style to "collim"  --  currently collimators have no slice number stored, so not too meaningful
    sliced_elem = create_thin_obj(thick_elem,1);
  }
  else if (strstr(thick_node->base_name,"solenoid"))    sliced_elem = create_thin_solenoid(thick_elem,1);    // create the first thin solenoid slice
  else if (strstr(thick_node->base_name,"elseparator")) sliced_elem = create_thin_elseparator(thick_elem,1); // create the first thin elseparator slice
  else // magnet which can be sliced to multipole
  {
    if(verbose>2) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << EOL;
    sliced_elem = create_sliced_magnet(thick_elem,1,ThickSLice); // get info from first slice
    if(ThickSLice) // create entry, body, exit pieces,  for bends or quadrupoles   --- if not yet existing
    {
      if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " ThickSLice, nslices=" << nslices << " create thick slices _en, _bo, _ex sliced_elem->name=" << sliced_elem->name << EOL;
      en =               create_thick_slice(thick_elem,0); // entry slice
      if(IsRbend && !MakeDipedge) add_half_angle_to(thick_elem,en,"e1");
      if(nslices>1) bo = create_thick_slice(thick_elem,1); // body slices,   last parameter = 1     since there will be only one type of body
      if(ThickSLice && IsQuad) ex=en; // for quad entry/exit are the same
      else ex =          create_thick_slice(thick_elem,nslices); // exit slice
      if(IsRbend && !MakeDipedge) add_half_angle_to(thick_elem,ex,"e2");
    }
  } // done with create_thin or create_thick.  Next is positioning slices as node

  if(verbose>1)
  {
    cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ <<  " thick " << thick_elem->name  << " of type " << thick_node->base_name << " done with element slice generation ";
    if(sliced_elem) cout << " sliced_elem->name=" << sliced_elem->name;
    if(en)          cout << " en->name=" << en->name; else cout << " en=" << en;
    if(bo)          cout << " bo->name=" << bo->name; else cout << " bo=" << bo;
    if(ex)          cout << " ex->name=" << ex->name; else cout << " ex=" << ex;
    cout << EOL;
  }

  if(verbose>2) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << EOL;
  command_parameter* length_param = return_param_recurse(const_cast<char*>("l"),thick_elem); // get original length, value or expression

  expression* at_expr = thick_node->at_expr;
  double at = thick_node->at_value;
  double length = thick_node->length; // direct curved thick_elem->length, for rbend longer than straight l

  if(verbose>2) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << EOL;
  double l_expr_val=length;
  expression* l_expr = NULL;
  if (length_param)
  {
    l_expr  = length_param->expr;
    if(l_expr) l_expr_val = expression_value(l_expr, 2);
  }

  if(verbose>2) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << EOL;
  double at_centre = get_node_pos(thick_node,thick_sequ); // for the moment just check - same with refer = center

  int middle=-1;
  if (nslices>1) middle = nslices/2; // used to determine after which slide to place the central marker

  if(verbose>1)
  {
    if(IsRbend &&ThickSLice) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " position sliced rbend, careful l_expr_val=" << l_expr_val << " " << my_dump_expression(l_expr) << EOL;
    cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " sliced_elem=" << left << setw(20) << sliced_elem->name << right
    << " thick_node->at_value=" << setw(8) << thick_node->at_value
    << " nslices=" << setw(2) << nslices << " at=" << setw(8) << at << " at_centre=" << setw(8) << at_centre << " length=" << setw(8) << length << " " << my_dump_element(sliced_elem) << " node:" << my_dump_node(thick_node);
  }

  if(IsBend)
  {
    if(EntryDipedge) // write start dipedge for dipoles --    for s - positions use the full curved length,  after rbend to sbend conversion, the sbend will be slightly longer and the start/end just outside the original rbend
      // make sure to write with enough precision, good choice is   set,format="20.14g";
    {
      if(UseDipedges) create_and_place_bend_node(thick_node,thick_node->at_expr,EntryDipedge,thin_sequ,-0.5,IsRbend); // subtract half of the length to be at start
    }
  }

  for (int i=0; i<nslices; ++i) // loop to place the nslices in the sequence
  {
    element *slice_i=NULL;
    if      (strstr(thick_node->base_name,"collimator"))  slice_i = create_thin_obj(thick_elem,i+1);
    else if (strstr(thick_node->base_name,"solenoid"))    slice_i = create_thin_solenoid(thick_elem,i+1);
    else if (strstr(thick_node->base_name,"elseparator")) slice_i = create_thin_elseparator(thick_elem,i+1);
    else // magnet which can be sliced to multipole
    {
      slice_i = create_sliced_magnet(thick_elem,i+1,ThickSLice); // create and place the multipole pieces
      if(ThickSLice) // fill space between slices
      {
        if(i==0) place_thick_slice(thick_elem,thick_node,thin_sequ,en,0,thin_style); // place entry  slice
        else     place_thick_slice(thick_elem,thick_node,thin_sequ,bo,i,thin_style); // place body/middle slice
        // place exit body after loop
      }
    }

    if(verbose>2) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " thick_node->base_name=" << thick_node->base_name << " i=" << i << " middle=" << middle << " nslices=" << nslices << " slice_i->name=" << slice_i->name << EOL;
    node* thin_node = new_elem_node(slice_i, thick_node->occ_cnt);
    thin_node->length   = 0.0;
    thin_node->from_name = buffer(thick_node->from_name);
    if (fabs(at_shift(nslices,i+1,local_thin_style))>0.0)
    {
      if (at_expr || l_expr)
      {
        thin_node->at_expr =
        compound_expr(at_expr,at,const_cast<char*>("+"),scale_expr(l_expr,at_shift(nslices,i+1,local_thin_style)),length*at_shift(nslices,i+1,local_thin_style));
      }
    }
    else
    {
      if (at_expr) thin_node->at_expr = clone_expression(at_expr);
    }
    thin_node->at_value = at + length*at_shift(nslices,i+1,local_thin_style);
    if (i==middle && !ThickSLice)
    {
      node* middle_marker_node=new_marker(thick_node,at,at_expr);
      add_node_at_end_of_sequence(middle_marker_node,thin_sequ);  // add a marker in the middle, except for thick nslices
      if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " thick_node " << thick_node->name << " place middle marker" << EOL;
    }

    if(verbose>1) { cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " nslices=" << nslices << " dump_node(thin_node):"; cout << my_dump_node(thin_node); }
    if(slice_i && !ThickSLice) add_node_at_end_of_sequence(thin_node,thin_sequ); // place thin multipole slices
  }

  if(ThickSLice)
  {
    place_thick_slice(thick_elem,thick_node,thin_sequ,ex,nslices,thin_style); // place exit slice
  }

  if(ExitDipedge && UseDipedges) if(UseDipedges) create_and_place_bend_node(thick_node,thick_node->at_expr,ExitDipedge,thin_sequ,0.5,IsRbend);  // write end dipedge for dipoles

} // SeqElList::slice_this_node()

element* SeqElList::create_thin_obj(const element* thick_elem, int slice_no) // creates the thin non-magnetic element - recursively, keeps original l as lrad
{
  element *thin_elem_parent = NULL;

  if (thick_elem == thick_elem->parent)
  {
	return NULL;
  }
  else
  {
	thin_elem_parent = create_thin_obj(thick_elem->parent,slice_no);
  }

  element *thin_elem;
  if( (thin_elem = theSliceList->find_slice(thick_elem,slice_no)) ) return thin_elem; // check to see if we've already done this one

  command* cmd = clone_command(thick_elem->def); // set up new multipole command
  const command_parameter* length_param = return_param_recurse(("l"),thick_elem);
  const int length_i = name_list_pos(const_cast<char*>("l"),thick_elem->def->par_names);
  const int lrad_i   = name_list_pos("lrad",thick_elem->def->par_names);
  if (length_param)
  {
	if (lrad_i > -1 && thick_elem->def->par_names->inform[lrad_i]>0)
	{ // already exists, replace lrad
	  cmd->par->parameters[lrad_i]->double_value = cmd->par->parameters[length_i]->double_value;
	  if (cmd->par->parameters[length_i]->expr)
	  {
		if (cmd->par->parameters[lrad_i]->expr)
		  delete_expression(cmd->par->parameters[lrad_i]->expr);
		cmd->par->parameters[lrad_i]->expr =
		clone_expression(cmd->par->parameters[length_i]->expr);
	  }
	}
	else // lrad does not yet exist
	{
	  if (name_list_pos("lrad",thick_elem->base_type->def->par_names) > -1)
	  { // add lrad only if allowed by element
		if (cmd->par->curr == cmd->par->max) grow_command_parameter_list(cmd->par);
		if (cmd->par_names->curr == cmd->par_names->max)
		  grow_name_list(cmd->par_names);
		cmd->par->parameters[cmd->par->curr] = clone_command_parameter(length_param);
		add_to_name_list(const_cast<char*>("lrad"),1,cmd->par_names);
		cmd->par->parameters[name_list_pos(const_cast<char*>("lrad"),cmd->par_names)]->expr =
		clone_expression(cmd->par->parameters[length_i]->expr);
		cmd->par->curr++;
	  }
	}
  }

  if (length_i > -1)
  {
	cmd->par->parameters[length_i]->double_value = 0;
	cmd->par->parameters[length_i]->expr = NULL;
  }
  int slices=1;
  if (strstr(thick_elem->base_type->name,"collimator")) slices = get_slices_from_elem(thick_elem);
  char* thin_name = NULL;
  if (slices==1 && slice_no==1) thin_name=buffer(CopToNewC_String(thick_elem->name));
  else thin_name=make_thin_name(thick_elem->name,slice_no);

  if (thin_elem_parent) thin_elem = make_element(thin_name,thin_elem_parent->name,cmd,-1);
  else  thin_elem = make_element(thin_name,thick_elem->base_type->name,cmd,-1);

  thin_elem->length = 0;
  thin_elem->bv = el_par_value(const_cast<char*>("bv"),thin_elem);

  theSliceList->put_slice(thick_elem,thin_elem);
  return thin_elem;
}

node* SeqElList::copy_thin(node* thick_node) // this copies an element node and sets the length to zero and radiation length to the length to be used for "copying" optically neutral elements
{
  node* thin_node = NULL;
  thin_node = clone_node(thick_node, 0);
  thin_node->length=0;
  thin_node->p_elem->length=0;
  if (el_par_value(const_cast<char*>("l"),thick_node->p_elem)>zero) thin_node->p_elem = create_thin_obj(thick_node->p_elem,1); // keep original length as lrad
  return thin_node;
}

void SeqElList::slice_node() // this decides how to split an individual node and sends it onto the thin_sequ builder
{
  node* thin_node;
  if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " thin_style=\"" << thin_style << "\"";
  if (thick_node->p_elem) // look at the element of this node to see what to do with slicing
  {
	if(verbose>1) cout << " now see what to do with thick_node=" << left << setw(20) << thick_node->name <<  " depending on its base=" << setw(20) << thick_node->base_name << right << EOL;
	const double eps=1.e-15; // used to check if a value is compatible with zero
	if ( fabs(el_par_value(const_cast<char*>("l"),thick_node->p_elem)) <eps ) // if the length is compatible with zero copy it directly
	{
	  if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " place copy of thick_node " << thick_node->name << EOL;
	  add_node_at_end_of_sequence(thin_node = copy_thin(thick_node),thin_sequ);
	}
	else if(strcmp(thick_node->base_name,"matrix") == 0) add_node_at_end_of_sequence(thick_node,thin_sequ); // Take matrix as it is, including any length
	else  // deal with elements which are just copied
	{
	  if (strcmp(thick_node->base_name,"marker") == 0      ||
		  strcmp(thick_node->base_name,"instrument") == 0  ||
		  strcmp(thick_node->base_name,"placeholder") == 0 ||
		  strcmp(thick_node->base_name,"hmonitor") == 0    ||
		  strcmp(thick_node->base_name,"vmonitor") == 0    ||
		  strcmp(thick_node->base_name,"monitor") == 0     ||
		  strcmp(thick_node->base_name,"vkicker") == 0     ||   // vkicker can have length, symplectic kick, leave as is
		  strcmp(thick_node->base_name,"hkicker") == 0     ||   // hkicker can have length, symplectic kick, leave as is
		  strcmp(thick_node->base_name,"kicker") == 0      ||
		  strcmp(thick_node->base_name,"tkicker") == 0     ||
		  strcmp(thick_node->base_name,"rfcavity") == 0    ||
		  strcmp(thick_node->base_name,"crabcavity") == 0
		  )
	  {
		if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " place copy of thick_node " << thick_node->name << EOL;
		add_node_at_end_of_sequence(thin_node = copy_thin(thick_node),thin_sequ);
		if (strcmp(thin_node->p_elem->base_type->name, "rfcavity") == 0 &&
			find_element(thin_node->p_elem->name, thin_sequ->cavities) == NULL)
		  add_to_el_list(&thin_node->p_elem, 0, thin_sequ->cavities, 0);     // special cavity
		if (strcmp(thin_node->p_elem->base_type->name, "crabcavity") == 0 &&
			find_element(thin_node->p_elem->name, thin_sequ->crabcavities) == NULL)
		  add_to_el_list(&thin_node->p_elem, 0, thin_sequ->crabcavities, 0); // special crab cavity
	  }
	  // now the element types that can be sliced
	  else if (strcmp(thick_node->base_name,"rbend") == 0       ||
			   strcmp(thick_node->base_name,"sbend") == 0       ||
			   strcmp(thick_node->base_name,"quadrupole") == 0  ||
			   strcmp(thick_node->base_name,"sextupole") == 0   ||
			   strcmp(thick_node->base_name,"octupole") == 0    ||
			   strcmp(thick_node->base_name,"solenoid") == 0    ||
			   strcmp(thick_node->base_name,"multipole") == 0   ||
			   strcmp(thick_node->base_name,"rcollimator") == 0 ||  // rcollimator has no slice  number, so will always be sliced to 1 slice
			   strcmp(thick_node->base_name,"ecollimator") == 0 ||  // ecollimator has no slice  number, so will always be sliced to 1 slice
			   strcmp(thick_node->base_name,"elseparator") == 0
			   )
	  {
		if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " " << thick_node->name <<  " base=" << setw(20) << thick_node->base_name << right << " slice_this_node" << EOL;
		slice_this_node(); // deal with elements and nodes that can be sliced, add the sliced element at_end_of_sequence
	  }
	  else if (strcmp(thick_node->base_name,"drift") == 0) ; // do nothing for drift
	  else // new elements not yet implemented for slicing - write a message and do a reasonable default action
	  {
		fprintf(prt_file, "Found unknown basename %s, place copy with length set to zero.\n",thick_node->base_name);
		add_node_at_end_of_sequence(copy_thin(thick_node),thin_sequ);
	  }
	}
  } // done with case where thick_node->p_elem  is defined
  else if (thick_node->p_sequ) // nested sequence (not flattened), slice and add the subsequence, this case is checked in fivecell.madx
  {
	sequence* sub_thin=slice_sequence(thin_style,thick_node->p_sequ);
	node* sub_node = new_sequ_node(sub_thin, thick_node->occ_cnt);
	sub_node->length = 0;
	sub_node->at_value = thick_node->at_value;
	if (sub_node->at_expr) sub_node->at_expr = clone_expression(thick_node->at_expr);
	if(verbose>1) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " place the sliced sub-sequence " << thick_node->p_sequ->name << EOL;
	add_node_at_end_of_sequence(sub_node,thin_sequ);
  }
  else fatal_error("node is not element or sequence",thick_node->base_name); // completely unknown, error
}

//--------  SequenceList
sequence* SequenceList::get_sequ(sequence* thick_sequ) // check if thick_sequ is already in my_sequ_list_vec
{
  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " my_sequ_list_vec.size()=" << my_sequ_list_vec.size() << EOL;
  for(unsigned int i=0; i<my_sequ_list_vec.size(); ++i)
  {
	if ( my_sequ_list_vec[i] == thick_sequ )
	{
	  if ( verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " found at i=" << i << EOL; // debug
	  return thick_sequ;
	}
  }
  return NULL;
}

void SequenceList::put_sequ(sequence* thick_sequ)
{
  my_sequ_list_vec.push_back(thick_sequ);
  if(verbose_fl()) cout << __FILE__<< " " << __FUNCTION__ << " line " << setw(4) << __LINE__ << " my_sequ_list_vec.size()=" << my_sequ_list_vec.size() << EOL;
  return;
}
