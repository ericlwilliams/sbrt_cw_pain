### Type1-Example1.py ###

from pynomo.nomographer import *

S_params={
        'u_min':230,
        'u_max':1,
        # 'function':lambda u:-(0.0174*u+log(0.0446)),
        'function':lambda u:-(0.0174*u),
        # 'title':r'$V_{99~\rm{Gy}_{2.1}}$',
         'title':r'$V_{30}$',
        'tick_levels':3,
        'tick_text_levels':1,
        'tick_side':'left',
                }

P_params={
        'u_min':0.9,
        'u_max':0.01,
        'function':lambda u:log(u/(1-u)),
        'title':r'$P$',
        'tick_levels':3,
        'tick_text_levels':1,
        'tick_side':'left',
                }

V_params={
        'u_min':48.3,
        'u_max':1,
        'function':lambda u:-(0.0407*u),
        'title':r'$V$',
        'tick_levels':3,
        'tick_text_levels':1,
                }

block_1_params={
        'block_type':'type_1',
        'width':10.0,
        'height':10.0,
        'f1_params':S_params,
        'f2_params':P_params,
        'f3_params':V_params,
               }

main_params={
        'filename':'test_cwp_1yr_nomo.pdf',
        'paper_height':10.0,
        'paper_width':10.0,
        'block_params':[block_1_params],
        'transformations':[('rotate',0.01),('scale paper',),
             ('polygon',)],
        'title_x':5.0,
        'title_y':-1.0,
        'title_box_width': 10.0,
        'title_str':r'$(S+0.64)^{0.58}(0.74V) = P$'
               }
Nomographer(main_params)
