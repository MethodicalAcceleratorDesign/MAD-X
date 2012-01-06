#ifndef MAD_GCST_H
#define MAD_GCST_H

/* Temporary file: global constants
   these constants will be split over their respective modules...
   and should be const pointers to constant values...
*/

extern const char *version_name;
extern const char *version_arch;
extern const char *version_type_dev;
extern const char *version_type_pro;
extern const char *version_date_mod;

extern char* const functs[];
extern const char op_string[];
extern char file_string[];        /* to avoid local in routine alias */
extern char vrai[];               /* to avoid local in routine alias */
extern char faux[];               /* to avoid local in routine alias */
extern const int n_match;         /* # of match token lists in cmd_match_base */
extern const int s_match[];       /* position of first token of command below */
extern const int t_match[];       /* order in which the commands are matched */
extern const char* cmd_match_base[];

/* aperture types and # of parameters, needed for twiss table */

extern char*  aperture_types[];

/*added 4, 3 and "racetrack" here */

extern int    aperture_npar[];

/* table descriptors: type 1 = int, type 2 = double, type 3 = string;
   internally, however, int are stored as double */

extern int    ap_table_types[];
extern char*  ap_table_cols[];
extern int    survey_table_types[];
extern char*  survey_table_cols[];
extern int    efield_table_types[];
extern char*  efield_table_cols[];
extern char*  sxf_table_names[];

extern int    twiss_opt_end;   /* last column filled by twiss module */
extern int    twiss_mult_end;  /* last multipole column filled by complete_twiss_table */
extern int    twiss_fill_end;  /* last standard column filled by complete_twiss_table */
/*== jln 11.11.2010 increased twiss_fill_end from 96 to 97 to accomodate for v_pos */

/* warning: modify routine complete_twiss_table in case of changes */
extern int    twiss_table_types[];
extern char*  twiss_table_cols[];
extern int    twiss_sector_table_types[];
extern char*  twiss_sector_table_cols[];
extern int    ptc_twiss_summary_table_types[];
extern char*  ptc_twiss_summary_table_cols[];

extern int    ibs_table_types[];
extern char*  ibs_table_cols[];
extern int    map_tab_types[];
extern char*  map_tab_cols[];
extern int    normal_res_types[];
extern char*  normal_res_cols[];

extern int    sodd_detune_5_types[];
extern char*  sodd_detune_5_cols[];
extern int    sodd_distort1_8_types[];
extern char*  sodd_distort1_8_cols[];
extern int    sodd_distort1_11_types[];
extern char*  sodd_distort1_11_cols[];
extern int    sodd_distort2_9_types[];
extern char*  sodd_distort2_9_cols[];

extern int    touschek_table_types[];
extern char*  touschek_table_cols[];
extern int    mon_table_types[];
extern char*  mon_table_cols[];
extern int    corr_table_types[];
extern char*  corr_table_cols[];

extern int    orbit_table_types[];
extern char*  orbit_table_cols[];
extern int    special_comm_cnt[];
extern char*  special_comm_desc[]; 
extern int    summ_table_types[]; 

extern char*  summ_table_cols[]; 

extern int    trackone_table_types[];
extern char*  trackone_table_cols[];
extern int    track_table_types[];
extern char*  track_table_cols[];
extern int    track_table_cols_len;
extern int    tracksumm_table_types[];
extern char*  tracksumm_table_cols[];
extern int    ptcnodetrack_table_types[];
extern char*  ptcnodetrack_table_cols[];


extern int    trackloss_table_types[];
extern char*  trackloss_table_cols[];

extern int    dynap_table_types[];
extern char*  dynap_table_cols[];
extern int    dynaptune_table_types[];
extern char*  dynaptune_table_cols[];

/* Definition of "select_ptc_normal" parameters for "ptc_normal"*/
extern char   names[PTC_NAMES_L][5];
extern char   atm[3][4];
extern char   atc[3][4];

extern char*  atc_type;
extern int    atc_flag;

extern char*  atm_type;
extern int    atm_flag;

#endif // MAD_GCST_H
