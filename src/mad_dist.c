#include "madx.h"
void setdistparameters(char *type, int cut_l, double* cuts, int index, int *dist_type, double *start_v, double *stop_v, double *coord_t){
	
	//1st First input (coordinate)
	//2nd coord type, 0 - action angle, 1 - normalized
	//3rd type 4 - uniform dist, 5 - normalized,  6-reyligh
	//5 start value for uniform, mean value for normalized and reyligh
	//6 sstop value for uniform, sigma for normalized and reyligh

	if(strcmp(type, "gauss")==0	){
		dist_type[index] = 6;
		dist_type[index+1] = 4;
		start_v[index] = 0;
		start_v[index+1] = 0;
		stop_v[index]  = 1;
		stop_v[index+1]  = twopi;
		for(int i=0; i < 6; i++)
			coord_t[i] = 0;

		if(cut_l==2){
			setactionanglecut(index, cuts[0], cuts[1]);	
		}
	}
	else if(strcmp(type, "uniform")==0){
		dist_type[index] = 4;
		dist_type[index+1] = 4;
		start_v[index] = pow(cuts[0],2);
		start_v[index+1] = 0;
		stop_v[index]  = pow(cuts[1],2);
		stop_v[index+1]  = twopi;
		for(int i=0; i < 6; i++)
			coord_t[i] = 0;
	}
	else if(strcmp(type, "fixed")==0){
		dist_type[index] = 0;
		dist_type[index+1] = 0;
		start_v[index] = 2; //2J
		start_v[index+1] = 0;
		stop_v[index]  = 0;
		stop_v[index+1]  = 0;
		for(int i=0; i < 6; i++)
			coord_t[i] = 0;
	}
	else if(strcmp(type, "zero")==0){
		dist_type[index] = 0;
		dist_type[index+1] = 0;
		start_v[index] = 0; //2J
		start_v[index+1] = 0;
		stop_v[index]  = 0;
		stop_v[index+1]  = 0;
		for(int i=0; i < 6; i++)
			coord_t[i] = 0;
	}
}


void pro_distribution(struct in_cmd* p){
	struct table* dist_t;
	char *table_name;
	char *type;
	double num ;
	double eigen [36];
	int dist_type [6];
	double start_v[6], stop_v[6], coord_t[6];
	int cut_l;
	double cuts[3];
	int const C_HOR = 0,C_VER = 2,C_LON = 4;
	double *xd, *pxd, *yd, *pyd, *td, *ptd;

	double npart_d = command_par_value("npart", p->clone);
	int npart = (int)npart_d;
	xd  = mymalloc("distribution", npart * sizeof *xd);
	pxd = mymalloc("distribution", npart * sizeof *xd);
	yd  = mymalloc("distribution", npart * sizeof *xd);
	pyd = mymalloc("distribution", npart * sizeof *xd);
	td  = mymalloc("distribution", npart * sizeof *xd);
	ptd = mymalloc("distribution", npart * sizeof *xd);
	


	//double px[npart], y[npart], py[npart], t[npart], pt[npart];
	
	initializedistribution(1);
	setemitt12(get_value("beam", "ex"), get_value("beam", "ey"));
	setemitt3(get_value("beam", "et"));
	sete0andmass0(get_value("beam", "energy"),get_value("beam", "mass") );
	settotalsteps(npart);
	if(command_par_value("use_intial", p->clone)!=0){

		double betx = command_par_value("betx", p->clone);
		double alfx = command_par_value("alfx", p->clone);
		double bety = command_par_value("bety", p->clone);
		double alfy = command_par_value("alfy", p->clone);
		settwisstas(betx, alfx, bety, alfy);

	}
	else{
		trbegn_(oneturnmat,eigen);
		settasmatrixtranspose(eigen);
		
	}

	type = command_par_string("horizontal", p->clone);
	cut_l = command_par_vector("cutsig_h", p->clone, cuts);
	setdistparameters(type, cut_l, cuts, C_HOR, dist_type, start_v,stop_v, coord_t);

	type = command_par_string("vertical", p->clone);
	cut_l = command_par_vector("cutsig_v", p->clone, cuts);
	setdistparameters(type, cut_l, cuts, C_VER, dist_type, start_v,stop_v, coord_t);

	type = command_par_string("longitudinal", p->clone);
	cut_l = command_par_vector("cutsig_l", p->clone, cuts);
	setdistparameters(type, cut_l, cuts, C_LON, dist_type, start_v,stop_v, coord_t);


	for(int i=0; i<6; i++){
		setscan_para_diagonal(i,coord_t[i],dist_type[i],start_v[i],stop_v[i]);
	}
	

	getunconvertedcoord(xd,pxd,yd,pyd,td,ptd, &npart);


	table_name = command_par_string("table", p->clone);
 	dist_t = make_table(table_name, "Distribution", dist_table_cols,dist_table_types, npart);
 	add_to_table_list(dist_t, table_register);


	for(int i=0; i< npart; i++){
		double_to_table_curr(table_name, "x", &xd[i]);
	 	double_to_table_curr(table_name, "px", &pxd[i]);
	 	double_to_table_curr(table_name, "y", &yd[i]);
	 	double_to_table_curr(table_name, "py", &pyd[i]);
	 	double_to_table_curr(table_name, "t", &td[i]);
	 	double_to_table_curr(table_name, "pt", &ptd[i]);
	 	num = (int)i;
	 	double_to_table_curr(table_name, "number", &num);
	 	augmentcountonly(table_name);	
	}

	  if (dist_t->header == NULL)  dist_t->header = new_char_p_array(40);

	//strncpy(tmp, t->org_sequ->name, NAME_L);
	//table_add_header(t, "@ SEQUENCE         %%%02ds \"%s\"", strlen(tmp),stoupper(tmp));
	char tmp[10];
	int i = get_string("beam", "particle", tmp);

	table_add_header(dist_t, "@ PARTICLE         %%%02ds \"%s\"", i, stoupper(tmp));
	table_add_header(dist_t, "@ MASS             %%le  %F", get_value("beam", "mass"));
	table_add_header(dist_t, "@ CHARGE           %%le  %F", get_value("beam", "charge"));
	table_add_header(dist_t, "@ ENERGY           %%le  %F", get_value("beam", "energy"));
	table_add_header(dist_t, "@ PC               %%le  %F", get_value("beam", "pc"));
	table_add_header(dist_t, "@ EX               %%le  %F", get_value("beam", "ex"));
	table_add_header(dist_t, "@ EY               %%le  %F", get_value("beam", "ey"));
	table_add_header(dist_t, "@ ET               %%le  %F", get_value("beam", "et"));
	type = command_par_string("horizontal", p->clone);
	table_add_header(dist_t, "@ DIST_TYPE_H         %%%02ds \"%s\"", strlen(type), stoupper(type));
	type = command_par_string("vertical", p->clone);
	table_add_header(dist_t, "@ DIST_TYPE_V         %%%02ds \"%s\"", strlen(type), stoupper(type));
	type = command_par_string("longitudinal", p->clone);
	table_add_header(dist_t, "@ DIST_TYPE_L         %%%02ds \"%s\"", strlen(type), stoupper(type));

	//free(xd);
	//free(pxd);
	//free(yd);
	//free(pyd);
	//free(td);
	//free(ptd);
	free_distribution(0); // free the memory in the distribution block

}