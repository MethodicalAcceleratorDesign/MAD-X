cl -c /Zm1000 madxn.c 
cl -c /Zm1000 gxx11psc.c 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix plot.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix gxx11ps.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix madxm.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix dynap.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix emit.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix twiss.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix match.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix matchsa.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix touschek.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix survey.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix trrun.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix util.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix orbf.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix ptc_dummy.F 
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix ibsdb.F 
lf95 -c a_scratch_size.f90           
lf95 -c b_da_arrays_all.f90          
lf95 -c c_dabnew.f90                 
lf95 -c d_lielib.f90                 
lf95 -c e_define_newda.f90           
lf95 -c f_newda.f90                  
lf95 -c g_newLielib.f90              
lf95 -c h_definition.f90             
lf95 -c i_tpsa.f90                   
lf95 -c j_tpsalie.f90                
lf95 -c k_tpsalie_analysis.f90       
lf95 -c l_complex_taylor.f90         
lf95 -c m_real_polymorph.f90         
lf95 -c n_complex_polymorph.f90      
lf95 -c o_tree_element.f90
lf95 -c Sa_extend_poly.f90           
lf95 -c Sb_1_pol_template.f90        
lf95 -c Sb_2_pol_template.f90        
lf95 -c Sc_euclidean.f90             
lf95 -c Sd_frame.f90                 
lf95 -c Se_status.f90                
lf95 -c Sf_def_all_kinds.f90         
lf95 -c Sg_0_fitted.f90              
lf95 -c Sg_1_template_my_kind.f90    
lf95 -c Sg_2_template_my_kind.f90    
lf95 -c Sh_def_kind.f90              
lf95 -c Si_def_element.f90           
lf95 -c Sj_elements.f90              
lf95 -c Sk_link_list.f90             
lf95 -c Sl_family.f90                
lf95 -c Sm_tracking.f90              
lf95 -c Sn_mad_like.f90              
lf95 -c So_fitting.f90               
lf95 -c Sp_keywords.f90              
lf95 -c madx_ptc_module.f90          
lf95 -c wrap.f90                     
lf95 -out madx madxm.obj madxn.obj dynap.obj emit.obj twiss.obj match.obj matchsa.obj touschek.obj survey.obj trrun.obj util.obj orbf.obj ibsdb.obj ptc_dummy.obj plot.obj gxx11ps.obj gxx11psc.obj 
lf95 -out madxdev madxm.obj madxn.obj dynap.obj emit.obj twiss.obj match.obj matchsa.obj touschek.obj survey.obj trrun.obj util.obj orbf.obj ibsdb.obj a_scratch_size.obj b_da_arrays_all.obj c_dabnew.obj d_lielib.obj e_define_newda.obj f_newda.obj g_newLielib.obj h_definition.obj i_tpsa.obj j_tpsalie.obj k_tpsalie_analysis.obj l_complex_taylor.obj madx_ptc_module.obj m_real_polymorph.obj n_complex_polymorph.obj o_tree_element.obj Sa_extend_poly.obj Sb_1_pol_template.obj Sb_2_pol_template.obj Sc_euclidean.obj Sd_frame.obj Se_status.obj Sf_def_all_kinds.obj Sg_0_fitted.obj Sg_1_template_my_kind.obj Sg_2_template_my_kind.obj Sh_def_kind.obj Si_def_element.obj Sj_elements.obj Sk_link_list.obj Sl_family.obj Sm_tracking.obj Sn_mad_like.obj So_fitting.obj Sp_keywords.obj wrap.obj plot.obj gxx11ps.obj gxx11psc.obj 
