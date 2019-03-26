#ifndef amp4hef_hpp
#define amp4hef_hpp

#define put_process      __amp4hef_MOD_put_process   
#define put_momenta      __amp4hef_MOD_put_momenta   
#define matrix_element_a __amp4hef_MOD_matrix_element_a
#define matrix_element_b __amp4hef_MOD_matrix_element_b
#define amplitude        __amp4hef_MOD_amplitude     

#include <complex>

extern "C" {

  void put_process( int *  ,int *  ,int *     ,int *   );
//                ( id     ,Ntotal ,Noffshell ,process ) 
//                ( output ,input  ,input     ,input   )

  void put_momenta( int * ,double * ,double *   );
//                ( id    ,momenta  ,directions ) 
//                ( input ,input    ,input      )

  void matrix_element_a( int * ,double *   );
  void matrix_element_b( int * ,double *   );
//                     ( id    ,ampSquared ) 
//                     ( input ,output     )

  void amplitude( int * ,std::complex<double> * ,int *    ,int *       );
//              ( id    ,ampValue               ,helicity ,permutation )
//              ( input ,output                 ,input    ,input       )

}

#endif
