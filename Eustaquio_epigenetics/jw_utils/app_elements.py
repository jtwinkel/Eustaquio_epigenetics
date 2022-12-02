from dash import dcc, html
import base64
from jw_utils import app_functions as afns
from jw_utils import file_utils as fu

def create_header(lab_name):
    image_filename = './assets/logo_simple.png'
    encoded_image = base64.b64encode(open(image_filename, 'rb').read())
    image = html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()),
                 height=100, width=100, id='logo',
                 style={'padding-top':0,
                        'padding-left':50,
                        'display':'inline-block',
                        })


    hd = html.Div([
            html.Div([
               image
            ], className='col-1',
                style={'width': '33%',
                    'display':'inline-block',
                    'float':'left',
                    }),
            html.Div([
                html.H1(['Disome project'], id='page-title')
            ], className="col-2",
                style={'width': '33%',
                    'display':'inline-block',
                    'vertical-align':'top',
                }
            ),
            html.Div([
                html.H1([lab_name], 
                    id='lab_name',
                    style={
                        'text-align':'right',
                        'padding-right':100,
                    }
                ),
            ], className="col-3",
                style={
                    'width': '33%',
                    'display':'inline-block',
                    'float':'right',
                    'vertical-align':'top',
                    'min-width':300
                }
            ),
        ], className="header-container",
            style={'width': '100%',
                   'padding-top':15,
                   'padding-bottom':0,
                   'overflow':'hidden',
                   'background-color':'rgb(248,248,248)',
                   'min-width':900})
    return hd
    
    
def create_feature_dropdown(path_to_gff):
    fd = html.Div([
        html.H4('Feature selector'),
        dcc.Dropdown(id = 'feature_dropdown',
            options = afns.add_selections_to_feature_dropdown(path_to_gff))
    ])
    return fd



def create_up_downstream_input():
    udi = html.Div([ 
        html.Div([
            html.H4('Upstream', style={'textAlign': 'left'}),
            dcc.Input(
                id="up_input",
                type="number",
                value=0,
            )
        ],  style = {
            'width':'40%',
            'display': 'inline-block'
            }
        ),
        html.Div([
            html.H4('Downstream', style={'textAlign': 'left'}),
            dcc.Input(
                id="down_input",
                type="number",
                value=0,
                )
        ],  style = {
            'width':'40%',
            'display': 'inline-block'
            }
        )
    ],  style = {
        'width':'100%',
        'overflow':'hidden',
        'min-width': 1500
        }
    )
    return udi
    
    
    
def create_experiment_selector(path_to_valcounts_dir):
    es = html.Div([
        html.H4('Experiment selector', style={'textAlign': 'left'}),
            dcc.Dropdown(id='experiment_selector',
                options = [{'label':afns.expname_from_path(path), 'value': afns.expname_from_path(path)} 
                                    for path in fu.get_filepaths_in_dir(path_to_valcounts_dir)],
                multi = True
            )
    ])
    return es
