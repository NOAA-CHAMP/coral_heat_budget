function stn = fixup_flkeys_mlrf1

  % E.g., http://tds.hycom.org/thredds/dodsC/flkeys.ascii?temperature[0:1463][0][244][298]

  stn = get_flkeys_hycom('mlrf1');

  stn.flkeys_hycom_seatemp.date = datenum(2008,1,1):(6/24):(datenum(2009,1,1)-(0.5/24));
  stn.flkeys_hycom_seatemp.data = [ 23.014206,23.06776,22.925781,23.123169,23.068495,23.05369,22.819693,22.624125,22.648241,22.65574,22.829426,22.963327,22.916216,23.002228,23.130863,23.14033,23.128942,23.098087,23.068857,23.05169,23.03335,23.026287,23.00099,23.35072,23.429491,24.023542,24.26767,24.31345,24.214827,23.98258,23.662542,23.74657,23.400414,23.266438,23.344492,23.423737,23.372469,23.359266,23.398155,23.399717,23.387814,23.356024,23.356747,23.341438,23.319818,23.280308,23.318205,23.338583,23.345747,23.391445,23.434141,23.425333,23.354233,23.36999,23.361788,23.184149,23.150589,23.087149,23.03177,22.956583,22.926392,22.966244,22.944565,22.9887,22.949144,22.977596,22.992767,23.129055,23.096539,23.073769,23.235695,23.307455,23.455439,23.911728,23.965534,23.961716,23.932201,23.893469,23.85573,23.722721,23.462328,23.368326,23.338856,23.268171,23.195234,23.154642,23.10661,23.074795,23.053036,23.017042,22.974176,22.970074,22.96,23.012985,23.07335,23.106668,23.143984,22.964994,22.873018,22.826757,22.790247,22.820498,22.780731,22.98171,23.05634,23.066175,22.968266,22.819202,22.690006,22.551235,22.489697,22.425652,22.426712,22.391016,22.428288,22.477701,22.409359,22.256138,22.40516,22.303968,22.321108,22.556078,22.776443,23.186253,23.191767,23.184858,23.223816,23.184456,23.043951,22.953722,22.657549,22.496502,22.76737,24.231655,23.996607,23.783674,23.537197,23.34112,23.540653,23.475801,23.466347,23.38858,23.323658,23.255892,23.227032,23.2993,23.217983,23.07526,23.085812,23.062336,22.902601,22.948881,23.191853,23.446991,23.527712,24.326237,23.975094,24.546633,24.152813,24.545761,24.54182,24.479656,24.297907,24.12582,23.96367,23.784607,23.727348,23.615307,23.514751,23.389904,23.359781,23.339275,23.222923,23.1893,23.224379,23.196917,22.973385,22.696344,22.776285,23.824942,23.808607,23.629208,23.98872,23.877947,23.620504,23.530243,23.443623,23.27485,23.250525,23.279196,23.19521,23.177017,23.235006,23.18305,23.011154,22.933838,23.08045,22.948605,22.967564,23.007822,22.90223,22.790312,22.748652,22.706608,22.660458,22.588816,22.56269,22.569515,22.614588,22.705883,22.906435,22.947126,22.955675,23.03274,23.075535,23.027046,23.076513,23.133495,23.33832,23.376604,23.561314,24.667255,23.997587,23.701769,23.937021,24.090105,24.414276,24.732521,24.690613,24.69848,24.557667,24.439175,23.764942,23.361668,23.00356,22.74334,22.571455,22.391356,22.209019,22.045847,21.983553,21.97286,22.048641,22.20341,22.165686,22.188063,22.215418,22.194788,22.226141,22.241621,22.221739,22.200415,22.18195,22.146782,22.177435,22.222689,22.255238,22.270824,22.380627,22.282518,22.463722,22.96922,23.372486,23.247204,23.364374,23.133368,22.533144,22.374336,22.324684,22.285248,22.250185,21.643412,21.846062,21.787619,21.677176,21.757423,21.703777,21.980696,22.307106,22.500307,22.468813,22.653727,22.607208,22.574972,22.552485,22.472351,22.562227,22.61776,22.427605,22.156496,22.245855,22.639385,22.774923,22.807327,22.918161,22.96221,22.950441,23.09355,22.573608,22.67351,22.909704,23.000248,23.222183,23.52399,23.669907,23.636269,23.119898,22.919409,22.730837,22.740528,22.904524,22.874043,23.151146,23.26406,23.26538,23.340664,23.343748,23.327646,23.280596,23.310514,23.383308,23.402533,23.35662,23.238058,23.110739,22.961637,23.029163,23.04738,23.060831,23.059565,23.0388,23.08581,23.160625,23.196934,23.211365,23.266804,23.197199,23.104395,22.909761,22.854649,22.797136,22.784643,22.774271,22.76359,22.75174,22.75656,22.760092,22.764421,22.776402,22.768394,22.761578,22.775059,22.79767,22.82842,22.859154,22.899158,23.025646,22.980457,23.018204,23.025791,23.029022,23.052425,23.06841,23.075914,23.076984,23.095299,23.147741,23.18788,23.244934,23.2383,23.150782,23.139986,23.142406,23.185722,23.118746,23.128662,23.156363,23.124817,23.29549,23.362299,23.342428,23.344439,23.269783,23.384642,23.43548,23.700226,24.033335,23.839607,23.966223,24.0766,24.299374,24.30492,24.179623,23.885386,24.247728,24.279116,24.354252,24.299343,24.227621,24.224665,24.18495,24.08229,24.119287,24.190063,24.194994,24.221224,24.321117,24.317167,24.315407,24.348936,24.400637,24.458925,24.50053,24.63975,24.690477,24.989477,24.985899,24.621574,24.426245,24.41161,24.477362,24.30506,24.153526,24.132254,24.14007,24.0184,23.900154,23.833294,23.711962,23.886436,23.977745,24.051674,24.038815,24.009958,24.016022,24.030426,24.067623,24.089014,24.030918,24.060503,24.012468,24.022686,23.961597,24.050104,24.082037,24.2393,24.248747,24.343948,24.311926,24.266811,24.134697,24.035763,24.081718,24.127794,24.129301,24.093676,24.202257,24.184868,24.18287,24.2106,24.169231,24.150812,24.218397,24.321081,24.634628,24.742527,24.844261,24.869087,24.834522,24.736334,24.805464,24.772123,24.78744,24.788599,24.740074,24.749123,24.88899,24.925558,24.99574,25.071371,25.125525,25.144798,25.115057,25.080261,24.974699,24.90172,24.870258,24.92121,24.921679,24.913881,24.925125,24.88135,24.82925,24.812298,24.934923,25.00403,25.050974,25.087524,25.088457,25.024145,25.079483,25.10385,25.166798,25.364737,25.549044,25.577394,25.257036,25.316359,25.265562,25.178722,25.251139,25.272781,25.161589,25.157434,25.205246,25.257149,25.401031,25.40781,25.465044,25.524637,25.682667,25.82268,26.007204,26.224966,26.366787,26.436903,26.485384,26.460152,26.679321,26.673927,26.41435,26.13372,26.106588,25.92955,25.818022,26.180712,26.21852,25.91846,25.989628,25.951183,26.025309,26.17864,26.174063,26.250774,26.190104,26.098864,26.005262,25.984411,26.069439,25.951216,25.954924,26.140038,26.28572,26.409988,26.583616,26.65219,26.678396,26.764458,26.779696,26.902899,27.042074,27.191753,27.145763,27.250664,27.37847,27.519323,27.853184,28.025064,28.018724,28.105944,27.974628,27.918375,27.89337,27.988806,27.616583,27.418552,27.412409,27.197458,27.233437,27.440979,27.696802,27.681707,27.868694,27.7486,27.49913,27.434855,27.34035,27.309757,27.227701,27.255177,27.2833,27.36674,27.361982,27.309643,27.328049,27.233433,27.16694,27.25704,27.527164,27.495024,27.390057,27.316536,27.28215,27.25693,27.206501,27.21606,27.181383,27.208307,27.205645,27.149725,27.104437,27.113382,27.249475,27.165228,27.078947,26.994867,27.063122,27.06915,27.000149,26.770535,26.863482,26.881598,26.9289,26.932913,26.93832,26.87142,26.836763,26.901232,27.051256,27.120253,27.104969,27.15644,27.1562,27.21829,27.476665,27.550882,27.50561,27.588682,27.594831,27.620058,27.56713,27.541294,27.532394,27.519863,27.526846,27.509262,27.480423,27.47608,27.531338,27.595135,27.63074,27.61262,27.614704,27.590525,27.6073,27.614819,27.696371,27.720558,27.688192,27.729692,27.912674,27.954689,27.974607,27.9457,27.924444,28.11796,28.105072,28.118328,28.148613,28.149452,28.124657,28.107775,28.01522,27.937939,27.900085,27.866665,27.85632,27.806183,27.784193,27.880194,27.811516,27.806795,27.822727,27.8145,27.867579,27.905523,27.928123,27.934715,27.907207,27.802082,28.023476,28.04983,28.12927,28.22529,28.241379,28.385496,28.599464,28.633837,28.660376,28.841663,29.044271,29.086397,28.92905,28.763723,28.714731,28.614912,28.392405,28.385925,28.230997,28.434414,28.28813,28.101233,28.092152,28.178259,28.15751,28.136074,28.28947,28.332497,28.341442,28.322994,28.36364,28.246197,28.233751,28.244776,28.51076,28.547844,28.459904,28.481972,28.567797,28.570719,28.52965,28.482538,28.26279,28.124004,27.795914,27.776741,28.29346,28.31285,28.106556,28.151093,28.206188,28.318777,28.360283,28.421734,28.363436,28.41123,28.48702,28.478586,28.462595,28.463644,28.34443,28.374779,28.342607,28.379663,28.354479,28.431505,28.41775,28.469408,28.33486,28.2683,28.329853,28.413862,28.357395,28.362745,28.39011,28.417953,28.391317,28.346712,28.317627,28.331518,28.268116,28.257612,28.343365,28.235659,28.288132,28.270727,28.380854,28.423166,28.418139,28.400946,28.405926,28.369146,28.27096,28.32988,28.266571,28.199158,28.152874,28.121286,28.113573,28.077456,28.046461,28.01508,28.189894,28.315634,28.360174,28.49537,28.596857,28.703995,28.737179,28.883217,28.879002,28.809757,28.796215,28.84633,28.847672,28.696022,28.533798,28.544441,28.550058,28.540966,28.485413,28.549942,28.618967,28.687857,28.728107,28.707155,28.706984,28.787626,28.771755,28.807629,28.854376,28.773308,28.850729,28.892239,28.942192,28.923342,28.863695,28.78751,28.716467,28.801107,28.867811,28.837872,28.80445,28.787388,28.90779,28.880722,28.814663,28.79442,28.69264,28.721592,28.724367,28.557245,28.490263,28.487402,28.472588,28.504139,28.5524,28.61422,28.77197,28.863976,28.906197,29.001999,29.180752,29.308355,29.397429,29.490093,29.589228,29.72344,29.71599,29.637701,29.588512,29.542313,29.395256,29.284132,29.272497,29.184042,29.079206,28.913586,29.064247,29.088211,29.122013,29.247526,29.246475,29.315338,29.333424,29.43372,29.510916,29.579779,29.675201,29.72514,29.764854,29.673033,29.764359,29.6849,29.636513,29.493553,29.536999,29.702251,29.643852,29.661686,29.694227,29.64294,29.6153,29.696577,29.717533,29.876509,29.973705,29.992916,30.00104,30.004747,30.008097,30.011456,30.102917,30.221867,30.241924,30.310493,30.34622,30.313103,30.300365,30.339603,30.500715,30.491209,30.46978,30.430296,30.431538,30.472115,30.41639,30.318707,30.247627,30.136126,29.996563,29.843206,29.787981,29.714645,29.539997,29.236801,29.276394,29.156584,28.770576,28.404097,28.802273,28.92004,28.90333,28.909111,28.910727,28.748272,28.658176,28.728437,28.6627,28.503716,28.506659,28.465023,28.562128,28.592443,28.553057,29.022596,29.110495,29.106712,29.109589,28.939648,28.83279,29.014765,29.459301,29.362125,29.280243,29.269827,29.344915,29.332775,29.423634,29.487776,29.58022,29.610466,29.606228,29.620672,29.668594,29.647137,29.633991,29.653027,29.582523,29.546215,29.471468,29.37917,29.292746,29.171726,29.085321,29.04601,29.087961,29.062359,29.066448,29.197006,29.209513,29.179785,29.202608,29.232342,29.175148,29.143177,29.13738,29.168377,29.06862,29.004902,29.098532,29.118095,29.063873,29.01644,29.055683,29.03224,29.0249,29.03014,29.028912,29.11777,29.159788,29.172216,29.227251,29.144709,29.085957,29.057415,29.054895,28.98855,28.916527,28.848675,28.725359,28.625399,28.503246,28.486563,28.358982,28.299608,28.214336,28.133043,28.052456,28.041775,28.027088,28.006813,28.043968,28.105154,28.090588,28.09916,28.161617,28.175323,28.169731,28.129759,28.100788,28.064411,28.02881,28.023073,27.9369,28.043833,28.248177,28.240046,28.23933,28.450405,28.584045,28.627008,28.71892,28.849361,28.936754,28.956455,28.850208,28.843866,28.872324,28.853388,28.851831,28.827168,28.736576,28.760527,28.7124,28.635384,28.606918,28.57546,28.524267,28.51224,28.522108,28.544626,28.486681,28.480917,28.45006,28.423151,28.35939,28.184185,28.260523,28.184547,28.189692,28.07638,28.251797,28.262365,28.20672,28.236256,28.19334,28.147804,28.126436,28.120636,28.132654,28.140137,28.153646,28.218317,28.184288,28.163792,28.169144,28.13345,28.062904,28.109295,28.121202,28.112408,28.136087,28.143585,28.155222,28.09507,28.060913,28.009085,27.9642,27.961542,27.927872,27.921837,27.85593,27.875414,27.968206,27.988476,27.930367,27.952614,27.965282,27.888186,27.815046,27.793873,27.746494,27.871536,28.099436,28.067791,28.052313,28.03339,28.045046,28.038967,28.022913,27.99967,27.977112,28.030828,28.03374,27.982086,27.856901,27.895582,27.954336,27.931398,27.874338,27.748632,27.795767,27.800467,27.775759,27.73559,27.841755,27.8639,27.90142,27.918055,27.925386,27.844,27.806665,27.80057,27.819338,27.8033,27.808016,27.78424,27.856318,27.826004,27.802711,27.812963,27.834063,27.664232,27.66144,27.703917,27.604221,27.653885,27.744621,27.718115,27.620651,27.61335,27.597153,27.57283,27.471571,27.51938,27.589073,27.646267,27.599453,27.523266,27.414253,27.491924,27.393974,27.35883,27.374622,27.373924,27.320948,27.267235,27.18873,27.152132,27.078968,27.059683,27.02437,27.015265,26.96508,26.94351,26.957933,26.938396,26.894901,26.82896,26.82163,26.774563,26.696217,26.680433,26.676928,26.646982,26.640411,26.700392,26.671167,26.617096,26.728745,26.69939,26.784567,26.747456,26.606537,26.383255,26.073263,25.851727,25.696362,25.538334,25.50561,25.551344,25.602144,25.512285,25.501856,25.622383,25.844929,26.071304,26.075024,25.996637,25.896172,25.83736,25.828377,25.776834,26.168941,26.135649,26.159193,25.881544,25.762127,25.747673,25.692612,25.722694,25.674475,25.645903,25.63901,25.640125,25.607458,25.559387,25.467596,25.407837,25.421026,25.407434,25.311714,25.261175,25.147572,25.283438,25.250847,25.287804,25.21563,25.279867,25.356796,25.359753,25.364323,25.343819,25.319077,25.317558,25.304506,25.285732,25.2542,25.23433,25.213402,25.18183,25.144178,25.159554,25.35472,25.395525,25.397272,25.358782,25.296602,25.19775,25.145092,25.070004,25.021622,25.02611,25.103533,25.1761,25.208967,25.260399,25.301632,25.276115,24.856764,24.740414,24.73209,24.86749,25.266134,25.108015,24.902992,24.76014,24.650936,24.521605,24.426416,24.282394,24.294876,24.298037,24.28327,24.156065,24.230097,24.315758,24.273424,24.10667,24.205357,24.145136,24.17063,24.120754,24.120361,24.042528,24.00518,23.927576,23.896269,23.848024,23.644817,23.739992,23.72016,23.677134,23.63528,23.601786,23.550468,23.493404,23.487104,23.414583,23.36899,23.390558,23.39186,23.344,23.296097,23.166794,23.114714,23.033224,22.940054,22.744207,22.676138,22.280128,21.781433,21.527958,21.43926,21.429594,21.421463,21.482769,21.42113,21.428417,21.321238,21.480019,21.247843,21.620157,21.554934,21.543312,21.554459,21.353048,21.46371,21.228748,21.633442,21.859781,22.053102,22.325539,22.795599,22.780125,22.371006,22.22286,22.521877,22.726637,22.75299,22.743372,22.77996,22.578743,22.17832,22.00606,21.944098,21.85803,21.794641,21.676855,21.511934,21.46353,21.92436,22.078289,21.98332,22.320093,22.22689,22.307434,22.462238,22.416695,22.454967,22.496584,22.46177,22.485502,22.356075,22.593967,22.791338,22.612892,22.550999,22.612497,22.663681,22.713558,22.50509,23.000067,22.906937,22.7332,22.68629,22.631992,22.643637,22.697136,22.812063,22.90092,23.15945,23.112644,22.962738,22.772419,22.713036,22.650175,22.603144,22.725788,22.808224,22.878138,22.90272,22.955074,22.853067,22.812338,22.87417,23.081491,23.225721,23.519726,23.710499,23.724495,23.751518,23.796282,23.704948,23.573769,23.442234,23.248962,22.960436,22.689558,22.830727,22.837212,22.792496,22.782831,22.728895,22.660103,22.547417,22.355743,22.302568,22.216227,22.253353,22.273045,22.268625,22.290203,22.334782,22.317543,22.472437,22.549257,22.648504,22.667467,22.831053,23.082426,23.333176,23.649656,23.949553,23.82028,23.70963,23.617989,23.559244,23.557676,23.468842,23.402018,23.327345,23.25928, ];

return;
