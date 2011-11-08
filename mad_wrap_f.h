#ifndef MAD_WRAP_F_H
#define MAD_WRAP_F_H

#ifndef WIN32

/* init fortran part */
#define mad_init_f            mad_init_f_

/* should work unchanged on _win32 using Lahey */
#define advance_node          advance_node_
#define advance_to_pos        advance_to_pos_
#define augment_count         augment_count_
#define augmentcountonly      augmentcountonly_
#define augmentcountmomtabs   augmentcountmomtabs_
#define char_from_table       char_from_table_  /* OB 2.4.2002 */
#define comment_to_table      comment_to_table_
#define comm_para             comm_para_
#define double_from_table     double_from_table_
#define string_from_table     string_from_table_ /* ETDA 8 nov 2004 */
#define double_to_table_row   double_to_table_row_ /* ETDA 11 nov 2004 */
#define double_table          double_table_    /* ETDA 25 aug 2004 */
#define double_to_table       double_to_table_

/* added by E. T. d'Amico on jan. 21st, 2004 */
#define interp_node           interp_node_
#define reset_interpolation   reset_interpolation_
#define embedded_twiss        embedded_twiss_
/* added by E. T. d'Amico on jan. 21stmay 19th, 2004 */
#define embedded_plot         embedded_plot_
/* end additions */

/*Piotr Skowronski (CERN)*/
#define gettrack gettrack_
#define deletetrackstrarpositions deletetrackstrarpositions_

#define copy_twiss_data       copy_twiss_data_
#define current_node_name     current_node_name_
#define element_name          element_name_
#define node_name             node_name_
#define frndm                 frndm_
#define el_par_vector         el_par_vector_
#define f_ctof                f_ctof_
#define get_disp0             get_disp0_
#define get_node_vector       get_node_vector_
#define get_option            get_option_
#define get_string            get_string_
#define get_title             get_title_
#define get_twiss_data        get_twiss_data_
#define get_variable          get_variable_
#define get_varstring         get_varstring_
#define get_vector            get_vector_
#define get_beam_value        get_beam_value_
#define get_value             get_value_
#define get_version           get_version_
#define grndm                 grndm_
#define headvalue             headvalue_
#define intrac                intrac_
#define mtcond                mtcond_
#define next_constraint       next_constraint_
#define next_constr_namepos   next_constr_namepos_
#define next_global           next_global_
#define getnumberoftracks     getnumberoftracks_
#define next_start            next_start_
#define next_vary             next_vary_
/* RDM 20.1.2006 BEGIN jacobian strategy (match) */
#define constraint_name       constraint_name_
#define vary_name             vary_name_
#define mtputconsname         mtputconsname_
/* RDM 20.1.2006 END jacobian strategy (match) */
#define node_al_errors        node_al_errors_
#define node_fd_errors        node_fd_errors_
#define node_string           node_string_
#define node_value            node_value_
#define plot_option           plot_option_
#define reset_count           reset_count_
#define restart_sequ          restart_sequ_
#define retreat_node          retreat_node_
#define sector_out            sector_out_
#define sequence_name         sequence_name_
#define set_option            set_option_
#define set_value             set_value_
#define set_variable          set_variable_
#define set_stringvar         set_stringvar_
#define spec_node_value       spec_node_value_
#define store_node_value      store_node_value_
#define store_node_vector     store_node_vector_
#define string_to_table       string_to_table_
#define table_length          table_length_
#define table_org             table_org_
#define table_range           table_range_
#define track_pteigen         track_pteigen_
#define vector_to_table       vector_to_table_
#define vdot                  vdot_
#define vmod                  vmod_
#define w_ptc_create_universe   w_ptc_create_universe_
#define w_ptc_create_layout     w_ptc_create_layout_
#define w_ptc_export_xml        w_ptc_export_xml_
#define w_ptc_move_to_layout    w_ptc_move_to_layout_
#define w_ptc_read_errors       w_ptc_read_errors_
#define w_ptc_refresh_k         w_ptc_refresh_k_
#define w_ptc_input             w_ptc_input_
#define w_ptc_align             w_ptc_align_
#define w_ptc_twiss             w_ptc_twiss_
#define w_ptc_normal            w_ptc_normal_
#define w_ptc_track             w_ptc_track_
#define w_ptc_start             w_ptc_start_
#define w_ptc_moments           w_ptc_moments_
#define w_ptc_initmoments       w_ptc_initmoments_
#define w_ptc_select            w_ptc_select_
#define w_ptc_writeparresults   w_ptc_writeparresults_
#define w_ptc_printframes       w_ptc_printframes_
#define w_ptc_printlayout_rootm w_ptc_printlayout_rootm_
#define w_ptc_eplacement        w_ptc_eplacement_
#define w_ptc_addknob           w_ptc_addknob_
#define w_ptc_addknob_i         w_ptc_addknob_i_
#define w_ptc_addmoment         w_ptc_addmoment_
#define w_ptc_getnmoments       w_ptc_getnmoments_
#define w_ptc_getmomentstabcol  w_ptc_getmomentstabcol_
#define w_ptc_setknobvalue      w_ptc_setknobvalue_
#define w_ptc_refreshtables     w_ptc_refreshtables_
#define w_ptc_getnfieldcomp     w_ptc_getnfieldcomp_
#define w_ptc_getsfieldcomp     w_ptc_getsfieldcomp_
#define w_ptc_setfieldcomp      w_ptc_setfieldcomp_
#define w_ptc_rviewer           w_ptc_rviewer_
#define w_ptc_script            w_ptc_script_
#define w_ptc_open_gino         w_ptc_open_gino_
#define w_ptc_addpush           w_ptc_addpush_
#define w_ptc_end               w_ptc_end_
#define w_ptc_dumpmaps          w_ptc_dumpmaps_
#define w_ptc_trackline         w_ptc_trackline_
#define w_ptc_track_everystep   w_ptc_track_everystep_
#define w_ptc_setdebuglevel     w_ptc_setdebuglevel_
#define w_ptc_enforce6d         w_ptc_enforce6d_
#define w_ptc_setaccel_method   w_ptc_setaccel_method_
#define w_ptc_setexactmis       w_ptc_setexactmis_
#define w_ptc_setradiation      w_ptc_setradiation_
#define w_ptc_setfringe         w_ptc_setfringe_
#define w_ptc_settotalpath      w_ptc_settotalpath_
#define w_ptc_settime           w_ptc_settime_
#define w_ptc_setnocavity       w_ptc_setnocavity_

#define seterrorflagfort        seterrorflagfort_  /*sets the dglobal error flag*/
#define geterrorflag            geterrorflag_  /*returns the dglobal error flag*/
#define getcurrentelementname   getcurrentelementname_
#define getcurrentcmdname       getcurrentcmdname_
#define stolower                stolower_
#define cf77flush               cf77flush_
#define select_ptc_idx          select_ptc_idx_  /* ETDA 10 nov 2004 */
#define result_from_normal      result_from_normal_ /* ETDA 11 nov 2004 */
#define make_map_table          make_map_table_ /* KZ 28.06.2005 table for maps */
#define minimum_acceptable_order minimum_acceptable_order_ /* ETDA 17 nov 2004 */
#define augmentfwarn            augmentfwarn_

#define makemomentstables  makemomentstables_

#define pro_input               pro_input_  /* AK 20070530*/

#define type_ofCall


#else /* WIN32 */

#define madx                  MADX

/* init fortran part */
#define mad_init_f            MAD_INIT_F

/* should work unchanged on _win32 using Lahey */
#define advance_node          ADVANCE_NODE
#define advance_to_pos        ADVANCE_TO_POS
#define augment_count         AUGMENT_COUNT
#define augmentcountonly      AUGMENTCOUNTONLY
#define augmentcountmomtabs   AUGMENTCOUNTMOMTABS
#define char_from_table       CHAR_FROM_TABLE  /* OB 2.4.2002 */
#define comment_to_table      COMMENT_TO_TABLE
#define copy_twiss_data       COPY_TWISS_DATA
#define current_node_name     CURRENT_NODE_NAME
#define comm_para             COMM_PARA
#define double_from_table     DOUBLE_FROM_TABLE
#define string_from_table     STRING_FROM_TABLE  /* ETDA 8 nov 2004 */
#define double_to_table_row   DOUBLE_TO_TABLE_ROW  /* ETDA 11 nov 2004 */
#define double_table          DOUBLE_TABLE    /* ETDA 25 aug 2004 */
#define double_to_table       DOUBLE_TO_TABLE

/* added by E. T. d'Amico on jan. 21st, 2004 */
#define interp_node           INTERP_NODE
#define reset_interpolation   RESET_INTERPOLATION
#define embedded_twiss        EMBEDDED_TWISS
/* added by E. T. d'Amico on jan. 21stmay 19th, 2004 */
#define embedded_plot         EMBEDDED_PLOT
/* end additions */

/*Piotr Skowronski (CERN)*/
#define gettrack              GETTRACK
#define deletetrackstrarpositions DELETETRACKSTRARPOSITIONS

#define element_name          ELEMENT_NAME
#define node_name             NODE_NAME
#define frndm                 FRNDM
#define el_par_vector         EL_PAR_VECTOR
#define f_ctof                F_CTOF
#define get_disp0             GET_DISP0
#define get_node_vector       GET_NODE_VECTOR
#define get_option            GET_OPTION
#define get_string            GET_STRING
#define get_title             GET_TITLE
#define get_twiss_data        GET_TWISS_DATA
#define get_variable          GET_VARIABLE
#define get_varstring         GET_VARSTRING
#define get_vector            GET_VECTOR
#define get_beam_value        GET_BEAM_VALUE
#define get_value             GET_VALUE
#define get_version           GET_VERSION
#define grndm                 GRNDM
#define headvalue             HEADVALUE
#define intrac                INTRAC
#define mtcond                MTCOND
#define next_constraint       NEXT_CONSTRAINT
#define next_constr_namepos   NEXT_CONSTR_NAMEPOS
#define next_global           NEXT_GLOBAL
#define getnumberoftracks     GETNUMBEROFTRACKS
#define next_start            NEXT_START
#define next_vary             NEXT_VARY
/* RDM 20.1.2006 BEGIN jacobian strategy (match) */
#define constraint_name       CONSTRAINT_NAME
#define vary_name             VARY_NAME
/* RDM 20.1.2006 END jacobian strategy (match) */
#define node_al_errors        NODE_AL_ERRORS
#define node_fd_errors        NODE_FD_ERRORS
#define node_string           NODE_STRING
#define node_value            NODE_VALUE
#define plot_option           PLOT_OPTION
#define reset_count           RESET_COUNT
#define restart_sequ          RESTART_SEQU
#define retreat_node          RETREAT_NODE
#define sector_out            SECTOR_OUT
#define sequence_name         SEQUENCE_NAME
#define set_option            SET_OPTION
#define set_value             SET_VALUE
#define set_variable          SET_VARIABLE
#define set_stringvar         SET_STRINGVAR
#define spec_node_value       SPEC_NODE_VALUE
#define store_node_value      STORE_NODE_VALUE
#define store_node_vector     STORE_NODE_VECTOR
#define string_to_table       STRING_TO_TABLE
#define table_length          TABLE_LENGTH
#define table_org             TABLE_ORG
#define table_range           TABLE_RANGE
#define track_pteigen         TRACK_PTEIGEN
#define vector_to_table       VECTOR_TO_TABLE
#define vdot                  VDOT
#define vmod                  VMOD
#define w_ptc_create_universe   W_PTC_CREATE_UNIVERSE
#define w_ptc_create_layout     W_PTC_CREATE_LAYOUT
#define w_ptc_export_xml        W_PTC_EXPORT_XML
#define w_ptc_move_to_layout    W_PTC_MOVE_TO_LAYOUT
#define w_ptc_read_errors       W_PTC_READ_ERRORS
#define w_ptc_refresh_k         W_PTC_REFRESH_K
#define w_ptc_input             W_PTC_INPUT
#define w_ptc_align             W_PTC_ALIGN
#define w_ptc_twiss             W_PTC_TWISS
#define w_ptc_normal            W_PTC_NORMAL
#define w_ptc_track             W_PTC_TRACK
#define w_ptc_start             W_PTC_START
#define w_ptc_moments           W_PTC_MOMENTS
#define w_ptc_initmoments       W_PTC_INITMOMENTS
#define w_ptc_select            W_PTC_SELECT
#define w_ptc_writeparresults   W_PTC_WRITEPARRESULTS
#define w_ptc_printframes       W_PTC_PRINTFRAMES
#define w_ptc_printlayout_rootm W_PTC_PRINTLAYOUT_ROOTM
#define w_ptc_eplacement        W_PTC_EPLACEMENT
#define w_ptc_addknob           W_PTC_ADDKNOB
#define w_ptc_addmoment         W_PTC_ADDMOMENT
#define w_ptc_getnmoments       W_PTC_GETNMOMENTS
#define w_ptc_getmomentstabcol  W_PTC_GETMOMENTSTABCOL
#define w_ptc_setknobvalue      W_PTC_SETKNOBVALUE
#define w_ptc_refreshtables     W_PTC_REFRESHTABLES
#define w_ptc_getnfieldcomp     W_PTC_GETNFIELDCOMP
#define w_ptc_getsfieldcomp     W_PTC_GETSFIELDCOMP
#define w_ptc_setfieldcomp      W_PTC_SETFIELDCOMP
#define w_ptc_rviewer           W_PTC_RVIEWER
#define w_ptc_script            W_PTC_SCRIPT
#define w_ptc_open_gino         W_PTC_OPEN_GINO
#define w_ptc_addpush           W_PTC_ADDPUSH
#define w_ptc_end               W_PTC_END
#define w_ptc_dumpmaps          W_PTC_DUMPMAPS
#define w_ptc_trackline         W_PTC_TRACKLINE
#define w_ptc_track_everystep   W_PTC_TRACK_EVERYSTEP
#define w_ptc_setdebuglevel     W_PTC_SETDEBUGLEVEL
#define w_ptc_enforce6d         W_PTC_ENFORCE6D
#define w_ptc_setaccel_method   W_PTC_SETACCEL_METHOD
#define w_ptc_setexactmis       W_PTC_SETEXACTMIS
#define w_ptc_setradiation      W_PTC_SETRADIATION
#define w_ptc_setfringe         W_PTC_SETFRINGE
#define w_ptc_settotalpath      W_PTC_SETTOTALPATH
#define w_ptc_settime           W_PTC_SETTIME
#define w_ptc_setnocavity       W_PTC_SETNOCAVITY

#define seterrorflagfort        SETERRORFLAGFORT  /*sets the dglobal error flag*/
#define geterrorflag            GETERRORFLAG  /*returns the dglobal error flag*/
#define getcurrentelementname   GETCURRENTELEMENTNAME
#define getcurrentcmdname       GETCURRENTCMDNAME

#define stolower                STOLOWER
#define cf77flush               CF77FLUSH
#define select_ptc_idx          SELECT_PTC_IDX  /* ETDA 10 nov 2004 */
#define result_from_normal      RESULT_FROM_NORMAL  /* ETDA 11 nov 2004 */
#define make_map_table          MAKE_MAP_TABLE  /* KZ 28.06.2005 table for maps */
#define minimum_acceptable_order MINIMUM_ACCEPTABLE_ORDER /* ETDA 17 nov 2004 */

#define augmentfwarn            AUGMENTFWARN

#define pro_input               PRO_INPUT   /* AK 20070530*/

#define type_ofCall _stdcall
#endif

#endif // MAD_WRAP_F_H

