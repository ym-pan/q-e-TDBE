  hP     k820309    -          22.0        BÌf                                                                                                          
       mbd_density.f90 MBD_DENSITY                                                     
                            @                              
       GEOM_T                      @                              
       DIAG                      @                              
       EIGVALSH INVERSE                                                       u #GET_DIAG_REAL    #GET_DIAG_COMPLEX    #MAKE_DIAG_REAL                                                          u #EIGVALSH_REAL    #EIGVALSH_COMPLEX 	                     @               Á           
     'x                   #COORDS    #LATTICE    #K_GRID    #CUSTOM_K_PTS    #PARALLEL_MODE    #GET_EIGS    #GET_MODES    #DO_RPA    #GET_RPA_ORDERS    #LOG    #TIMER )   #EXC ;   #FREQ @   #GAMM D   #REAL_SPACE_CUTOFF E   #REC_SPACE_CUTOFF F   #PARAM G   #IDX R   #INIT W   #DESTROY Z   #SIZ ]   #HAS_EXC `   #CLOCK c              $                                                           
            &                   &                                                       $                                         `                 
            &                   &                                                       $                                          À                             &                                                       $                                                         
            &                   &                                                        $                                  
       h                                                                                          Cauto                                       $                                   t                                                                                        @                   $                                   x                                                                                        B                   $                                   |                                                                                        D                   $                                         	                                                                                  F                    $                                                
       #LOGGER_T                   @  @                              '                    #LEVEL    #PRINTER    #INFO    #DEBUG    #WARN !   #ERROR %                $                                                                                                                               1    #         @   $                     d                           #STR                                                    y                                                         
                                                  1 1         À    $                                              #LOGGER_INFO    #         @     @                                                #THIS    #STR              
                                                    #LOGGER_T              
                                                    1 1         À    $                                              #LOGGER_DEBUG    #         @     @                                                #THIS    #STR               
                                                    #LOGGER_T              
                                                     1 1         À    $                            !                  #LOGGER_WARN "   #         @     @                            "                    #THIS #   #STR $             
                                 #                   #LOGGER_T              
                                $                    1 1         À    $                            %                  #LOGGER_ERROR &   #         @     @                            &                    #THIS '   #STR (             
                                 '                   #LOGGER_T              
                                (                    1               $                              )     à                    #CLOCK_T *                  @  @              E           *     'à                    #ACTIVE +   #LEVEL ,   #TIMESTAMPS -   #COUNTS .   #LEVELS /   #INIT 0   #CLOCK 4   #PRINT 8                $                              +                                                                                    ÿÿÿÿÿÿÿÿ   Q                   $                              ,                                                                                                0                $                             -                                         &                                                      $                             .            P                             &                                                       $                              /                                         &                                           1         À    $                            0                  #CLOCK_INIT 1   #         @     @                            1                    #THIS 2   #N 3             
                                2     à               #CLOCK_T *             
                                  3           1         À    $                            4                  #CLOCK_CLOCK 5   #         @     @                            5                    #THIS 6   #ID 7             
                                6     à               #CLOCK_T *             
                                  7           1         À    $                            8                  #CLOCK_PRINT 9   #         @     @                            9                    #THIS :             
                                :     à               #CLOCK_T *                 $                              ;     Ì       x             #EXCEPTION_T <                  @  @                         <     'Ì                    #CODE =   #ORIGIN >   #MSG ?                $                              =                                                                                                 0                 $                             >     2                                                                             3       °              C(unknown)                                                                          $                             ?            6                                                                             ²              C                                                                                                                                                                                      $                              @            H                   #QUAD_PT_T A             &                                                          @  @                         A     '                    #VAL B   #WEIGHT C                 $                             B                
                 $                             C               
                $                             D             
                                                 
                                 0D0                  $                             E              
                 $                             F               
                 $                              G     8       ¨             #PARAM_T H                     @                         H     '8              	      #DIPOLE_CUTOFF I   #EWALD_REAL_CUTOFF_SCALING J   #EWALD_REC_CUTOFF_SCALING K   #K_GRID_SHIFT L   #EWALD_ON M   #ZERO_NEGATIVE_EIGVALS N   #RPA_RESCALE_EIGS O   #RPA_ORDER_MAX P   #N_FREQ Q                $                             I               
                                                   
                  vú¤@                         $                             J              
                                                 
                       ð?        1D0                 $                             K              
                                                 
                       ð?        1D0                 $                             L              
                                                  
                       à?        0.5D0                 $                              M                                                                                    ÿÿÿÿÿÿÿÿ   %                   $                              N     $                                                                                         &                   $                              O     (                                                                                         '                   $                              P     ,                                                                           
               10                 $                              Q     0       	                                                                                    15                  $                              R            à             #ATOM_INDEX_T S                  @  @              D           S     '                    #I_ATOM T   #J_ATOM U   #N_ATOMS V               $                              T                                          &                                                       $                              U            H                             &                                                         $                              V                  1         À    $                            W                  #GEOM_INIT X   #         @     @                            X                    #THIS Y             
                                Y     x              #GEOM_T 
   1         À    $                            Z                  #GEOM_DESTROY [   #         @     @                            [                    #THIS \             
                                \     x              #GEOM_T 
   1         À    $                          ]                  #GEOM_SIZ ^   %         @   @                          ^                           #THIS _             
                                 _     x             #GEOM_T 
   1         À    $                           `                  #GEOM_HAS_EXC a   %         @   @                           a                           #THIS b             
                                 b     x             #GEOM_T 
   1         À    $                            c                  #GEOM_CLOCK d   #         @     @                            d                    #THIS e   #ID f             
                                e     x              #GEOM_T 
             
                                  f           (        `                               g                    /                
    #A h   #EXC i     p        H r j     7
S
O
 p        j                      j                        n                                       1            p          H r j     7
S
O
 p        j                      j                        n                                          1              H r j     7
S
O
 p        j                      j                        n                                      2                H r j     7
S
O
 p        j                      j                        n                                          1              H r j     7
S
O
 p        j                      j                        n                                      2                                                        
                                h                   
 -             &                   &                                                                                      i     Ì               #EXCEPTION_T <                                                k                                                                                                     l     
                   
                  -DTû!	@        (        `   @                                               	                
    #A m   p          H r n     7
S
O
 p        j                      j                        n                                    1                H r n     7
S
O
 p        j                      j                        n                                          1                                                
                                m                   
              &                   &                                           (        `   @                                                                   #A o   p          H r p     7SO p        j                      j                        n                                    1                H r p     7SO p        j                      j                        n                                          1                                                
                                o                    
             &                   &                                           (        `   @                                                              
    #D q     p        H r r     7
S
O
 p        j            j                                  p          H r r     7
S
O
 p        j            j                                    H r r     7
S
O
 p        j            j                                      H r r     7
S
O
 p        j            j                                    H r r     7
S
O
 p        j            j                                                              
                                q                   
              &                                           (        `   @                                              4                
    #A s   #EXC t   #DESTROY u   p          H r v     7
S
O
 p        j                      j                        n                                      1                H r v     7
S
O
 p        j                      j                        n                                          1                                              
                                s                   
 0             &                   &                                                                                      t     Ì               #EXCEPTION_T <             
                                 u           (        `   @                           	                    9                
    #A w   #EXC x   #DESTROY y   p          H r z     7SO p        j                      j                        n                                      1                H r z     7SO p        j                      j                        n                                          1                                              
                                w                    5             &                   &                                                                                      x     Ì               #EXCEPTION_T <             
                                 y           (        `                                {                    	                
    #GEOM |   #PTS }   #CHARGES ~   #MASSES    #OMEGAS    p          H r      7
S
O
 p        j                      j                        n                                           2                H r      7
S
O
 p        j                      j                        n                                          2                                         
                                  |     x             #GEOM_T 
          0  
 @                              }                   
              &                   &                                                     
                                 ~                   
              &                                                     
                                                    
              &                                                     
                                                    
              &                                           (        `                                                                    
    #GEOM    #PTS    #CHARGES    #MASSES    #OMEGAS    #MODES    p          H r      7
S
O
 p        j                      j                        n                                        2                H r      7
S
O
 p        j                      j                        n                                          2                                            
                                       x             #GEOM_T 
          0  
 @                                                 
 
             &                   &                                                     
                                                    
              &                                                     
                                                    
              &                                                     
  @                                                 
              &                                                     
  @                                                 
              &                   &                                                                                          SIZE                                                SIZE               @                           z     SIZE               @                           v     SIZE               @                           r     SIZE               @                           p     SIZE               @                           n     SIZE               @                           j     SIZE        $      fn#fn    Ä   @   J   MBD_CONSTANTS      G   J  MBD_GEOM    K  E   J  MBD_LINALG      Q   J  MBD_LAPACK $   á  }       gen@DIAG+MBD_LINALG (   ^  i       gen@EIGVALSH+MBD_LAPACK     Ç        GEOM_T+MBD_GEOM '   H  ¬   a   GEOM_T%COORDS+MBD_GEOM (   ô  ¬   a   GEOM_T%LATTICE+MBD_GEOM '         a   GEOM_T%K_GRID+MBD_GEOM -   4  ¬   a   GEOM_T%CUSTOM_K_PTS+MBD_GEOM .   à  Ç   a   GEOM_T%PARALLEL_MODE+MBD_GEOM )   §  ¤   a   GEOM_T%GET_EIGS+MBD_GEOM *   K  ¤   a   GEOM_T%GET_MODES+MBD_GEOM '   ï  ¤   a   GEOM_T%DO_RPA+MBD_GEOM /   	  ¤   a   GEOM_T%GET_RPA_ORDERS+MBD_GEOM $   7
  ^   a   GEOM_T%LOG+MBD_GEOM #   
        LOGGER_T+MBD_UTILS )   '  ¥   a   LOGGER_T%LEVEL+MBD_UTILS 3   Ì  ½      LOGGER_T%PRINTER+MBD_UTILS=PRINTER /     L   a   LOGGER_T%PRINTER%STR+MBD_UTILS (   Õ  Y   a   LOGGER_T%INFO+MBD_UTILS &   .  [      LOGGER_INFO+MBD_UTILS +     V   a   LOGGER_INFO%THIS+MBD_UTILS *   ß  L   a   LOGGER_INFO%STR+MBD_UTILS )   +  Z   a   LOGGER_T%DEBUG+MBD_UTILS '     [      LOGGER_DEBUG+MBD_UTILS ,   à  V   a   LOGGER_DEBUG%THIS+MBD_UTILS +   6  L   a   LOGGER_DEBUG%STR+MBD_UTILS (     Y   a   LOGGER_T%WARN+MBD_UTILS &   Û  [      LOGGER_WARN+MBD_UTILS +   6  V   a   LOGGER_WARN%THIS+MBD_UTILS *     L   a   LOGGER_WARN%STR+MBD_UTILS )   Ø  Z   a   LOGGER_T%ERROR+MBD_UTILS '   2  [      LOGGER_ERROR+MBD_UTILS ,     V   a   LOGGER_ERROR%THIS+MBD_UTILS +   ã  L   a   LOGGER_ERROR%STR+MBD_UTILS &   /  ]   a   GEOM_T%TIMER+MBD_GEOM "     ¯      CLOCK_T+MBD_UTILS )   ;  ¤   a   CLOCK_T%ACTIVE+MBD_UTILS (   ß  ¥   a   CLOCK_T%LEVEL+MBD_UTILS -        a   CLOCK_T%TIMESTAMPS+MBD_UTILS )        a   CLOCK_T%COUNTS+MBD_UTILS )   ¬     a   CLOCK_T%LEVELS+MBD_UTILS '   @  X   a   CLOCK_T%INIT+MBD_UTILS %     Y      CLOCK_INIT+MBD_UTILS *   ñ  U   a   CLOCK_INIT%THIS+MBD_UTILS '   F  @   a   CLOCK_INIT%N+MBD_UTILS (     Y   a   CLOCK_T%CLOCK+MBD_UTILS &   ß  Z      CLOCK_CLOCK+MBD_UTILS +   9  U   a   CLOCK_CLOCK%THIS+MBD_UTILS )     @   a   CLOCK_CLOCK%ID+MBD_UTILS (   Î  Y   a   CLOCK_T%PRINT+MBD_UTILS &   '  R      CLOCK_PRINT+MBD_UTILS +   y  U   a   CLOCK_PRINT%THIS+MBD_UTILS $   Î  a   a   GEOM_T%EXC+MBD_GEOM &   /  o      EXCEPTION_T+MBD_UTILS +     ¥   a   EXCEPTION_T%CODE+MBD_UTILS -   C  ï   a   EXCEPTION_T%ORIGIN+MBD_UTILS *   2  S  a   EXCEPTION_T%MSG+MBD_UTILS %     £   a   GEOM_T%FREQ+MBD_GEOM $   (  e      QUAD_PT_T+MBD_UTILS (     H   a   QUAD_PT_T%VAL+MBD_UTILS +   Õ  H   a   QUAD_PT_T%WEIGHT+MBD_UTILS %     §   a   GEOM_T%GAMM+MBD_GEOM 2   Ä  H   a   GEOM_T%REAL_SPACE_CUTOFF+MBD_GEOM 1      H   a   GEOM_T%REC_SPACE_CUTOFF+MBD_GEOM &   T   ]   a   GEOM_T%PARAM+MBD_GEOM !   ±         PARAM_T+MBD_GEOM /   Á!  ¤   a   PARAM_T%DIPOLE_CUTOFF+MBD_GEOM ;   e"  §   a   PARAM_T%EWALD_REAL_CUTOFF_SCALING+MBD_GEOM :   #  §   a   PARAM_T%EWALD_REC_CUTOFF_SCALING+MBD_GEOM .   ³#  ©   a   PARAM_T%K_GRID_SHIFT+MBD_GEOM *   \$  ¤   a   PARAM_T%EWALD_ON+MBD_GEOM 7    %  ¤   a   PARAM_T%ZERO_NEGATIVE_EIGVALS+MBD_GEOM 2   ¤%  ¤   a   PARAM_T%RPA_RESCALE_EIGS+MBD_GEOM /   H&  ¦   a   PARAM_T%RPA_ORDER_MAX+MBD_GEOM (   î&  ¦   a   PARAM_T%N_FREQ+MBD_GEOM $   '  b   a   GEOM_T%IDX+MBD_GEOM '   ö'  u      ATOM_INDEX_T+MBD_UTILS .   k(     a   ATOM_INDEX_T%I_ATOM+MBD_UTILS .   ÿ(     a   ATOM_INDEX_T%J_ATOM+MBD_UTILS /   )  H   a   ATOM_INDEX_T%N_ATOMS+MBD_UTILS %   Û)  W   a   GEOM_T%INIT+MBD_GEOM #   2*  R      GEOM_INIT+MBD_GEOM (   *  T   a   GEOM_INIT%THIS+MBD_GEOM (   Ø*  Z   a   GEOM_T%DESTROY+MBD_GEOM &   2+  R      GEOM_DESTROY+MBD_GEOM +   +  T   a   GEOM_DESTROY%THIS+MBD_GEOM $   Ø+  V   a   GEOM_T%SIZ+MBD_GEOM "   .,  Z      GEOM_SIZ+MBD_GEOM '   ,  T   a   GEOM_SIZ%THIS+MBD_GEOM (   Ü,  Z   a   GEOM_T%HAS_EXC+MBD_GEOM &   6-  Z      GEOM_HAS_EXC+MBD_GEOM +   -  T   a   GEOM_HAS_EXC%THIS+MBD_GEOM &   ä-  X   a   GEOM_T%CLOCK+MBD_GEOM $   <.  Z      GEOM_CLOCK+MBD_GEOM )   .  T   a   GEOM_CLOCK%THIS+MBD_GEOM '   ê.  @   a   GEOM_CLOCK%ID+MBD_GEOM #   */  -      INVERSE+MBD_LAPACK %   W3  ¤   a   INVERSE%A+MBD_LAPACK '   û3  Y   a   INVERSE%EXC+MBD_LAPACK !   T4  p       DP+MBD_CONSTANTS !   Ä4  p       PI+MBD_CONSTANTS )   45  õ     GET_DIAG_REAL+MBD_LINALG +   )7  ¤   a   GET_DIAG_REAL%A+MBD_LINALG ,   Í7  õ     GET_DIAG_COMPLEX+MBD_LINALG .   Â9  ¤   a   GET_DIAG_COMPLEX%A+MBD_LINALG *   f:  Ë     MAKE_DIAG_REAL+MBD_LINALG ,   1=     a   MAKE_DIAG_REAL%D+MBD_LINALG )   ½=       EIGVALSH_REAL+MBD_LAPACK +   È?  ¤   a   EIGVALSH_REAL%A+MBD_LAPACK -   l@  Y   a   EIGVALSH_REAL%EXC+MBD_LAPACK 1   Å@  @   a   EIGVALSH_REAL%DESTROY+MBD_LAPACK ,   A       EIGVALSH_COMPLEX+MBD_LAPACK .   C  ¤   a   EIGVALSH_COMPLEX%A+MBD_LAPACK 0   ´C  Y   a   EIGVALSH_COMPLEX%EXC+MBD_LAPACK 4   D  @   a   EIGVALSH_COMPLEX%DESTROY+MBD_LAPACK (   MD  &      EVAL_MBD_NONINT_DENSITY -   sF  T   a   EVAL_MBD_NONINT_DENSITY%GEOM ,   ÇF  ¤   a   EVAL_MBD_NONINT_DENSITY%PTS 0   kG     a   EVAL_MBD_NONINT_DENSITY%CHARGES /   ÷G     a   EVAL_MBD_NONINT_DENSITY%MASSES /   H     a   EVAL_MBD_NONINT_DENSITY%OMEGAS %   I  1      EVAL_MBD_INT_DENSITY *   @K  T   a   EVAL_MBD_INT_DENSITY%GEOM )   K  ¤   a   EVAL_MBD_INT_DENSITY%PTS -   8L     a   EVAL_MBD_INT_DENSITY%CHARGES ,   ÄL     a   EVAL_MBD_INT_DENSITY%MASSES ,   PM     a   EVAL_MBD_INT_DENSITY%OMEGAS +   ÜM  ¤   a   EVAL_MBD_INT_DENSITY%MODES *   N  =      EVAL_MBD_INT_DENSITY%SIZE -   ½N  =      EVAL_MBD_NONINT_DENSITY%SIZE 6   úN  =      EIGVALSH_COMPLEX%SIZE+MBD_LAPACK=SIZE 3   7O  =      EIGVALSH_REAL%SIZE+MBD_LAPACK=SIZE 4   tO  =      MAKE_DIAG_REAL%SIZE+MBD_LINALG=SIZE 6   ±O  =      GET_DIAG_COMPLEX%SIZE+MBD_LINALG=SIZE 3   îO  =      GET_DIAG_REAL%SIZE+MBD_LINALG=SIZE -   +P  =      INVERSE%SIZE+MBD_LAPACK=SIZE 