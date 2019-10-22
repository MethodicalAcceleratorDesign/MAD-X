#ifndef MAD_PTC_H
#define MAD_PTC_H

// types

struct in_cmd;

// interface

int   select_ptc_idx(void);
void  select_ptc_normal(struct in_cmd*);
void  pro_ptc_normal(void);
void  pro_ptc_twiss(void);
int   minimum_acceptable_order(void);
void  pro_ptc_read_errors(void);
void  pro_ptc_refresh_k(void);
void  pro_ptc_create_layout(void);
void  ptc_track_end(void);
void  ptc_track_observe(struct in_cmd*);
void  ptc_dumpmaps(struct in_cmd*);
void  ptc_oneturnmap(struct in_cmd*); // empty...
void  pro_ptc_track(struct in_cmd*);
void  pro_ptc_open_gino(struct in_cmd*);
void  pro_ptc_script(struct in_cmd*);
int   makemomentstables(void);
int   pro_ptc_select_moment(struct in_cmd*);
int   pro_ptc_moments(struct in_cmd*);
void  pro_ptc_select(struct in_cmd*);
void  pro_ptc_setfieldcomp(struct in_cmd*);
void  pro_ptc_setknobvalue(struct in_cmd*);
void  pro_ptc_knob(struct in_cmd*);
void  pro_ptc_varyknob(struct in_cmd*);
void  pro_ptc_eplacement(struct in_cmd*);
void  pro_ptc_export_xml(struct in_cmd*);
void  pro_ptc_printframes(struct in_cmd*);
void  pro_ptc_printparametric(struct in_cmd*);
void  pro_ptc_setswitch(struct in_cmd*);
void  pro_ptc_enforce6d(struct in_cmd*);
void  pro_ptc_trackline(struct in_cmd*);
void  ptc_putbeambeam(struct in_cmd*);


void  augmentcountmomtabs(double* s);

#endif // MAD_PTC_H

