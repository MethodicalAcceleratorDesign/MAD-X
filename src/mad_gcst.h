#ifndef MAD_GCST_H
#define MAD_GCST_H

/* Temporary file: global constants
   these constants will be split over their respective modules...
   and should be const pointers to constant values...
*/

extern const char* const version_name;
extern const char* const version_arch;
extern const char* const version_ostype;
extern const char* const version_date;

extern const char* const functs[];
extern const char* const op_string;
extern const char* const file_string;        /* to avoid local in routine alias */
extern const char* const vrai;               /* to avoid local in routine alias */
extern const char* const faux;               /* to avoid local in routine alias */
extern const int n_match;         /* # of match token lists in cmd_match_base */
extern const int s_match[];       /* position of first token of command below */
extern const int t_match[];       /* order in which the commands are matched */
extern const char* const cmd_match_base[];

/* aperture types and # of parameters, needed for twiss table */

extern const char* const aperture_types[];

/*added 4, 3 and "racetrack" here */

extern const int         aperture_npar[];

/* table descriptors: type 1 = int, type 2 = double, type 3 = string;
   internally, however, int are stored as double */

extern const int         ap_table_types[];
extern const char* const ap_table_cols[];
extern const int         survey_table_types[];
extern const char* const survey_table_cols[];
extern const int         efield_table_types[];
extern const char* const efield_table_cols[];
extern const char* const sxf_table_names[];

extern const int         twiss_opt_end;   /* last column filled by twiss module */
extern const int         twiss_mult_end;  /* last multipole column filled by complete_twiss_table */
extern const int         twiss_fill_end;  /* last standard column filled by complete_twiss_table */
/*== jln 11.11.2010 increased twiss_fill_end from 96 to 97 to accomodate for v_pos */

/* warning: modify routine complete_twiss_table in case of changes */
extern const int         twiss_table_types[];
extern const char* const twiss_table_cols[];
extern const int         twiss_sector_table_types[];
extern const char* const twiss_sector_table_cols[];
extern const int         ptc_twiss_summary_table_types[];
extern const char* const ptc_twiss_summary_table_cols[];

extern const int         bb6d_ixy_types[];
extern const char* const bb6d_ixy_cols[];

extern const int         ibs_table_types[];
extern const char* const ibs_table_cols[];
extern const int         map_tab_types[];
extern const char* const map_tab_cols[];
extern const int         normal_res_types[];
extern const char* const normal_res_cols[];

extern const int         sodd_detune_5_types[];
extern const char* const sodd_detune_5_cols[];
extern const int         sodd_distort1_8_types[];
extern const char* const sodd_distort1_8_cols[];
extern const int         sodd_distort1_11_types[];
extern const char* const sodd_distort1_11_cols[];
extern const int         sodd_distort2_9_types[];
extern const char* const sodd_distort2_9_cols[];

extern const int         touschek_table_types[];
extern const char* const touschek_table_cols[];
extern const int         mon_table_types[];
extern const char* const mon_table_cols[];
extern const int         corr_table_types[];
extern const char* const corr_table_cols[];

extern const int         orbit_table_types[];
extern const char* const orbit_table_cols[];
extern const int         special_comm_cnt[];
extern const char* const special_comm_desc[]; 
extern const int         summ_table_types[]; 

extern const char* const summ_table_cols[]; 

extern const int         trackone_table_types[];
extern const char* const trackone_table_cols[];
extern const int         track_table_types[];
extern const char* const track_table_cols[];
extern const int         track_table_cols_len;
extern const int         tracksumm_table_types[];
extern const char* const tracksumm_table_cols[];
extern const int         ptcnodetrack_table_types[];
extern const char* const ptcnodetrack_table_cols[];

extern const int         trackloss_table_types[];
extern const char* const trackloss_table_cols[];

extern const char* const dist_table_cols[];
extern const int         dist_table_types[];

extern const int         dynap_table_types[];
extern const char* const dynap_table_cols[];
extern const int         dynaptune_table_types[];
extern const char* const dynaptune_table_cols[];

/* Definition of "select_ptc_normal" parameters for "ptc_normal"*/
extern const char* const names[];
extern const char* const atm[];
extern const char* const atc[];

extern const char*  atc_type;
extern       int    atc_flag;

extern const char*  atm_type;
extern       int    atm_flag;

extern const char* const nonlin_table_cols[];  
extern const int    nonlin_table_types[];

#endif // MAD_GCST_H
