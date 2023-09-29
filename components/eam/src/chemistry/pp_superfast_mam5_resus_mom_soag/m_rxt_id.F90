










                                                                                
      module m_rxt_id                                                           
                                                                                
      implicit none                                                             
                                                                                
      integer, parameter :: rid_jo1dU =    1                                    
      integer, parameter :: rid_jo2_b =    2                                    
      integer, parameter :: rid_jh2o2 =    3                                    
      integer, parameter :: rid_jch2o_a =    4                                  
      integer, parameter :: rid_jch2o_b =    5                                  
      integer, parameter :: rid_jch3ooh =    6                                  
      integer, parameter :: rid_jc2h5ooh =    7                                 
      integer, parameter :: rid_jno2 =    8                                     
      integer, parameter :: rid_jno3_a =    9                                   
      integer, parameter :: rid_jno3_b =   10                                   
      integer, parameter :: rid_jn2o5_a =   11                                  
      integer, parameter :: rid_jn2o5_b =   12                                  
      integer, parameter :: rid_jhno3 =   13                                    
      integer, parameter :: rid_jho2no2_a =   14                                
      integer, parameter :: rid_jho2no2_b =   15                                
      integer, parameter :: rid_jch3cho =   16                                  
      integer, parameter :: rid_jpan =   17                                     
      integer, parameter :: rid_jacet =   18                                    
      integer, parameter :: rid_jmvk =   19                                     
      integer, parameter :: rid_jsoa_a1 =   20                                  
      integer, parameter :: rid_jsoa_a2 =   21                                  
      integer, parameter :: rid_jsoa_a3 =   22                                  
      integer, parameter :: rid_uci1 =   23                                     
      integer, parameter :: rid_uci2 =   24                                     
      integer, parameter :: rid_uci3 =   25                                     
      integer, parameter :: rid_lco_h =   26                                    
      integer, parameter :: rid_lco_ho2 =   27                                  
      integer, parameter :: rid_lh2_ho2 =   28                                  
      integer, parameter :: rid_lch4 =   29                                     
      integer, parameter :: rid_lc2h6 =   30                                    
      integer, parameter :: rid_lc3h8 =   31                                    
      integer, parameter :: rid_lc2h4_oh =   32                                 
      integer, parameter :: rid_lc2h4_o3 =   33                                 
      integer, parameter :: rid_lisop_o3 =   34                                 
      integer, parameter :: rid_lisop_oh =   35                                 
      integer, parameter :: rid_lch2o =   36                                    
      integer, parameter :: rid_lo3_oh =   37                                   
      integer, parameter :: rid_po3_oh =   38                                   
      integer, parameter :: rid_lo3_ho2 =   39                                  
      integer, parameter :: rid_lho2_oh =   40                                  
      integer, parameter :: rid_uci4 =   41                                     
      integer, parameter :: rid_uci5 =   42                                     
      integer, parameter :: rid_ph2o2 =   43                                    
      integer, parameter :: rid_lh2o2 =   44                                    
      integer, parameter :: rid_lo3_no =   45                                   
      integer, parameter :: rid_lno_ho2 =   46                                  
      integer, parameter :: rid_lo3_no2 =   47                                  
      integer, parameter :: rid_lno3_oh =   48                                  
      integer, parameter :: rid_lno3_no =   49                                  
      integer, parameter :: rid_lhno4 =   50                                    
      integer, parameter :: rid_lhno3 =   51                                    
      integer, parameter :: rid_uci6 =   52                                     
      integer, parameter :: rid_lno2_oh =   53                                  
      integer, parameter :: rid_HO2NO2f =   54                                  
      integer, parameter :: rid_N2O5f =   55                                    
      integer, parameter :: rid_PANf =   56                                     
      integer, parameter :: rid_uci7 =   57                                     
      integer, parameter :: rid_uci8 =   58                                     
      integer, parameter :: rid_uci9 =   59                                     
      integer, parameter :: rid_lch3o2_ho2 =   60                               
      integer, parameter :: rid_lch3o2_no =   61                                
      integer, parameter :: rid_lch3o2 =   62                                   
      integer, parameter :: rid_lch3ooh =   63                                  
      integer, parameter :: rid_lc2h5o2_no =   64                               
      integer, parameter :: rid_lc2h5o2 =   65                                  
      integer, parameter :: rid_lc2h5o2_ch3 =   66                              
      integer, parameter :: rid_lc2h5o2_ho2 =   67                              
      integer, parameter :: rid_lc2h5ooh_a =   68                               
      integer, parameter :: rid_lc2h5ooh_b =   69                               
      integer, parameter :: rid_lch3cho_oh =   70                               
      integer, parameter :: rid_lch3cho_no3 =   71                              
      integer, parameter :: rid_lch3co3_no =   72                               
      integer, parameter :: rid_lch3co3_ch3 =   73                              
      integer, parameter :: rid_lch3co3 =   74                                  
      integer, parameter :: rid_lch3coch3_a =   75                              
      integer, parameter :: rid_lch3coch3_b =   76                              
      integer, parameter :: rid_lroho2_no =   77                                
      integer, parameter :: rid_lroho2_ho2 =   78                               
      integer, parameter :: rid_lroho2_ch3o2 =   79                             
      integer, parameter :: rid_lisopo2_no =   80                               
      integer, parameter :: rid_lisopo2_ho2 =   81                              
      integer, parameter :: rid_lisopo2_ch3 =   82                              
      integer, parameter :: rid_lmvkmacr_o3 =   83                              
      integer, parameter :: rid_lmvkmacr_oh =   84                              
      integer, parameter :: rid_lmvko2_no =   85                                
      integer, parameter :: rid_lmvko2_ho2 =   86                               
      integer, parameter :: rid_usr_e90 =   87                                  
      integer, parameter :: rid_ldms_oh =   88                                  
      integer, parameter :: rid_usr_DMS_OH =   89                               
      integer, parameter :: rid_usr_SO2_OH =   90                               
      integer, parameter :: rid_ldms_no3 =   91                                 
      integer, parameter :: rid_vsoag0_oh =   92                                
      integer, parameter :: rid_vsoag15_oh =   93                               
      integer, parameter :: rid_vsoag24_oh =   94                               
      integer, parameter :: rid_visop_oh =   95                                 
      integer, parameter :: rid_vc10h16_oh =   96                               
      integer, parameter :: rid_visop_o3 =   97                                 
      integer, parameter :: rid_vc10h16_o3 =   98                               
      integer, parameter :: rid_visop_no3 =   99                                
      integer, parameter :: rid_vc10h16_no3 =  100                              
      integer, parameter :: rid_vsoag35_oh =  101                               
      integer, parameter :: rid_vsoag34_oh =  102                               
      integer, parameter :: rid_vsoag33_oh =  103                               
      integer, parameter :: rid_vsoag32_oh =  104                               
                                                                                
                                                                                
      end module m_rxt_id                                                       
