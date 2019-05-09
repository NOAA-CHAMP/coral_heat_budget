function stn = fixup_flkeys_sanf1

  % E.g., http://tds.hycom.org/thredds/dodsC/flkeys.ascii?temperature[0:1463][0][184][148]

  stn = get_flkeys_hycom('sanf1');

  stn.flkeys_hycom_seatemp.date = datenum(2008,1,1):(6/24):(datenum(2009,1,1)-(0.5/24));
  stn.flkeys_hycom_seatemp.data = [ 24.961746,24.945852,24.86733,24.899147,24.81016,24.655815,24.379004,24.264883,24.071028,23.69102,23.443134,23.302309,23.229792,23.22626,23.137543,23.07337,23.221045,23.19191,23.265987,23.172413,22.844564,23.298065,23.207571,23.37152,23.276646,23.344229,23.343596,23.320292,23.280804,23.215342,23.423119,23.440119,23.33326,23.185335,23.412743,23.436502,23.424713,23.399757,23.421843,23.427107,23.42199,23.42019,23.507206,23.528254,23.590582,23.565899,23.52148,23.518671,23.561962,23.588564,23.634329,23.665949,23.653137,23.66061,23.629383,23.560793,23.37898,23.301744,23.245255,23.208672,23.167171,23.158108,23.205734,23.24523,23.255785,23.303055,23.330719,23.340275,23.391016,23.44165,23.4549,23.517998,23.562426,23.581343,23.495695,23.397547,23.376682,23.3411,23.256786,23.06892,23.07888,23.105923,23.135393,23.129494,23.1534,23.177288,23.173353,23.17019,23.15493,23.1444,23.128296,23.076311,22.986214,22.976072,23.00903,22.64219,21.97369,21.9882,22.173298,22.334656,22.504845,22.459557,22.570225,22.718357,22.864872,22.950775,22.948431,22.773512,21.847519,21.277945,21.89054,22.55241,22.46919,22.47307,22.43848,22.420963,22.41305,22.440899,22.436363,22.382442,22.459045,22.440681,22.225773,22.021646,21.869633,22.088287,22.119408,22.255957,22.270788,22.254976,21.313316,21.625416,21.964474,21.897394,21.737867,21.85315,22.485565,22.5002,22.729498,22.787622,22.476841,22.452017,22.80646,22.8736,22.811934,23.089071,23.249758,23.214933,23.509254,23.61422,23.722855,23.746807,23.773542,23.738363,22.306782,23.442461,23.707108,23.677963,23.589178,23.072826,22.260546,22.161077,22.133724,22.26944,22.291069,22.502983,22.821878,22.927732,22.899775,22.78878,22.812815,22.985525,23.074986,23.05264,23.213165,23.228134,23.39428,23.289701,23.182058,22.931444,23.004242,23.144104,23.128225,23.121902,23.132256,23.185577,23.192385,23.203081,23.253633,23.270887,23.269806,23.210693,23.182716,23.340818,23.360128,23.399717,23.536879,23.550194,23.47079,23.10705,22.849188,22.870743,22.875385,22.81079,22.845333,22.84835,22.903368,22.959068,23.031067,23.11306,23.25568,23.48342,23.69974,23.89483,23.68475,23.896473,23.756151,23.854311,23.87435,23.905075,23.733402,23.773777,23.93885,23.462784,23.53824,23.448107,23.480072,23.424576,23.281586,23.10506,22.722712,22.404667,22.36301,22.01799,21.912933,21.777693,21.47435,21.385052,21.33445,21.368366,21.583416,21.794266,22.06247,22.038427,21.821877,21.769602,21.741514,21.686914,21.813437,21.860504,22.07755,22.831957,22.981554,22.772882,22.821526,22.569725,22.655285,22.579103,22.548897,22.446852,22.418158,22.486824,22.535547,22.686909,22.623777,22.591967,22.792849,23.072338,23.04229,22.76323,22.596754,22.491838,22.497366,22.264608,22.011335,21.876644,21.661364,21.595985,21.56006,21.422997,21.333113,21.54455,21.762085,22.815754,23.075495,22.928165,22.564991,22.299402,21.956717,21.547962,21.491758,21.653885,21.628485,21.6605,21.710901,21.75802,21.852154,22.001759,22.072191,22.038317,22.448555,22.323833,22.24021,22.379711,22.23255,22.350073,22.082516,22.18615,22.31179,22.906443,23.429098,23.383179,23.346695,23.37829,23.397354,23.3529,23.396196,23.395021,23.454117,23.423483,23.4347,23.542252,23.469912,23.428411,23.435478,23.491028,23.52727,23.625298,23.697386,23.61965,23.649092,23.753296,23.866957,23.717808,23.604849,23.586775,23.414293,23.218576,22.749592,22.524527,22.425552,22.3632,22.39235,22.503557,22.618008,22.715693,22.892694,23.003483,22.950203,22.937935,22.943829,23.0331,22.96657,23.008608,22.965233,22.985449,23.066954,23.12778,23.141218,23.173273,23.133457,23.146078,23.23854,23.269491,23.309782,23.34881,23.368858,23.397781,23.421906,23.453562,23.478487,23.488916,23.493883,23.473587,23.47079,23.454979,23.589432,23.741629,23.985683,24.173527,24.209656,24.198421,24.277006,24.327223,24.38942,24.397104,24.35665,24.349281,24.358887,24.384007,24.399263,24.475489,24.49851,24.463066,24.478228,24.425682,24.46701,24.396915,24.352985,24.339561,24.33283,24.307486,24.31046,24.298843,24.29351,24.28145,24.273619,24.29519,24.349874,24.397142,24.466505,24.54331,24.797777,24.949118,25.05078,25.118637,24.827936,24.676958,24.785158,25.03339,25.292238,24.87732,24.486324,24.466457,24.595106,24.537952,24.226465,24.11412,24.066864,24.00506,23.958809,23.927937,23.866947,23.80512,23.75968,23.7642,23.786282,23.80065,23.782135,23.816439,23.796764,23.813915,23.994833,24.17541,24.375467,24.469587,24.414473,24.45123,24.637976,24.441376,24.189249,24.343592,24.57509,24.536509,24.346745,24.413593,24.495111,24.490833,24.313438,24.377655,24.278162,24.161575,24.138472,24.152617,24.17853,24.231552,24.27287,24.359053,24.588963,24.542578,24.58952,24.625786,24.630466,24.615591,24.650408,24.695736,24.736664,24.871143,24.964605,24.858881,24.841553,24.886929,25.00001,25.053661,25.31935,25.333664,24.98772,24.840662,24.8079,24.846823,24.870907,24.890463,24.880285,24.873682,24.87478,24.940483,24.889359,24.878242,24.910704,24.965725,25.065681,25.208853,25.275331,25.467926,25.795856,25.622627,25.434858,25.947939,26.161236,26.273735,26.28061,26.485704,26.365597,26.352324,26.364214,26.10659,26.451666,26.331242,25.904316,26.259594,26.486671,26.49575,26.827305,27.929327,28.337267,28.376802,28.386408,27.88406,28.045576,27.638617,27.680069,27.396263,26.921045,27.074936,27.859867,27.640085,27.39647,27.172216,26.959131,26.819126,26.745045,26.746927,26.717531,26.661297,26.627876,26.55615,26.409729,26.388247,26.721025,26.8257,27.060425,27.279547,27.142729,27.038322,27.201527,26.456373,26.830168,27.21658,27.253746,27.383825,27.565453,27.535042,27.12465,26.815132,26.734953,27.092873,27.053236,26.877544,27.133741,27.496033,27.289694,27.323082,27.52353,27.759293,27.517162,27.221203,27.595898,27.239492,27.196032,27.26592,27.612251,27.586802,27.62311,27.55614,27.325773,27.322886,27.4157,27.223907,27.061235,26.84988,26.805641,26.859186,27.01501,27.184229,27.171183,27.138996,27.073648,27.0208,26.978846,26.930191,26.912771,26.91512,26.92527,26.936884,26.918102,26.97884,26.971508,26.925783,27.02093,27.083288,27.128614,27.202143,27.340038,27.32861,27.331264,27.328674,27.357054,27.529268,27.82484,28.061216,28.199656,28.126669,28.070545,28.128927,28.456852,28.370386,27.956215,27.963024,27.792519,27.687817,27.654589,27.639946,27.651686,27.67679,27.690035,27.7157,27.700584,27.663424,27.665432,27.678852,27.696945,27.698566,27.587347,27.471083,27.502193,27.5301,27.625393,27.913036,28.427742,28.586254,28.581104,28.358934,28.556543,28.708715,29.021523,29.110296,29.060167,29.093618,29.123034,29.098202,29.099266,28.980091,28.831604,28.797743,28.760199,28.585869,28.67411,28.677431,28.551884,28.372818,28.372818,28.371433,28.393015,28.432158,28.304989,28.196892,28.263119,28.520958,27.986359,27.799274,27.808409,27.980562,27.586948,27.33812,27.630966,27.891706,27.716583,27.55346,27.679403,27.699339,27.865725,28.0085,28.098854,27.710762,28.202972,28.59571,28.720284,28.54901,28.68043,28.801113,28.795643,28.156181,28.599867,28.614853,28.53783,28.400623,28.451597,28.38032,28.378056,28.35048,28.302816,28.452406,28.507278,28.478128,28.40005,28.394732,28.355364,28.319761,28.482504,28.684492,28.709167,28.654818,28.581837,28.571474,28.53533,28.504614,28.444624,28.491808,28.50746,28.493155,28.429588,28.510038,28.530132,28.477163,28.494465,28.52072,28.56873,28.578535,28.597036,28.624111,28.6671,28.664486,28.655136,28.637007,28.624826,28.592545,28.555815,28.563267,28.507746,28.533327,28.489887,28.486359,28.481789,28.485126,28.484228,28.501085,28.510258,28.518442,28.535772,28.534815,28.515446,28.536324,28.506458,28.467098,28.434566,28.41968,28.40929,28.421957,28.426254,28.426338,28.447353,28.451145,28.454823,28.457966,28.466597,28.467781,28.484812,28.505129,28.570793,28.66257,28.787844,28.728626,28.518576,28.477213,28.386513,28.274414,28.26467,28.331617,28.337872,28.432198,28.41447,28.44393,28.555721,28.735956,28.885925,29.198206,29.331366,29.43424,29.295553,28.999178,28.875408,28.887156,28.881739,28.888512,28.936108,28.974596,28.914509,28.991488,28.987305,29.011564,29.044855,29.009054,29.002728,29.02392,29.017721,28.987701,28.918518,28.977776,28.925856,28.921692,28.88877,28.927114,29.017912,28.964968,28.921555,28.923626,28.943707,28.926605,28.92263,28.867573,28.870407,28.891499,28.859615,28.841845,28.987078,28.973732,28.66337,28.657097,28.760338,28.77965,28.82035,28.845167,28.920776,28.909452,28.967146,28.999233,29.152216,29.162489,29.177574,29.178423,29.225645,29.341768,29.329355,29.447533,29.456482,29.447279,29.467878,29.492926,29.466005,29.45738,29.56176,29.53564,29.485228,29.492155,29.549686,29.509064,29.527538,29.57537,29.55898,29.579418,29.540281,29.506304,29.420956,29.403658,29.415388,29.456125,29.502644,29.527067,29.594757,29.709377,29.874731,29.98488,30.027254,30.073273,30.131098,30.255684,30.220385,30.305515,30.367058,30.375715,30.477,30.220997,30.1712,30.206911,30.34166,30.266727,30.11202,30.258642,30.396511,30.298855,30.18926,30.299416,30.509928,30.522692,30.489042,30.572523,30.45307,30.128216,30.180233,30.472036,30.4872,30.512367,30.568192,30.648878,30.531492,30.641695,30.564358,30.241762,30.13164,30.0242,29.892654,29.855406,29.743471,29.576782,29.438354,29.420502,29.54113,29.125408,28.764534,28.956903,29.154219,29.118187,28.92732,29.132423,29.264616,29.260645,29.382757,29.496809,29.504265,29.420454,29.452879,29.625721,29.775764,29.854914,29.978287,29.887545,29.923437,29.914873,29.960438,30.105104,30.139347,30.128416,30.241638,30.187784,30.056297,29.782764,29.709976,29.757677,29.755898,29.839813,29.857313,29.721025,29.634718,29.75896,30.029507,30.120016,30.054178,29.971067,29.8464,29.799799,29.704653,29.62812,29.57762,29.529226,29.472624,29.42652,29.402704,29.382565,29.437313,29.646856,29.722939,29.51518,29.43631,29.447563,29.418926,29.345488,29.229118,29.142313,29.104574,29.080576,29.05738,29.022503,29.01901,29.06499,29.084097,29.234545,29.263435,29.403294,29.147846,28.943436,28.9226,28.98238,29.002565,29.239292,29.408525,29.365276,29.315622,29.213354,29.079737,28.931366,28.7101,28.588749,28.4739,28.352858,28.299093,28.328579,28.42326,28.46529,28.490686,28.489935,28.486277,28.49306,28.481026,28.472713,28.389078,28.32166,28.399954,28.376781,28.407026,28.377401,28.418362,28.475073,28.812449,28.934158,28.999353,28.861025,28.693941,28.779411,28.75606,29.178783,29.602598,29.280378,29.029758,28.897802,28.629478,28.843676,29.195446,28.95281,28.825718,28.917067,28.865261,28.88602,28.8075,28.782858,28.729433,28.729921,28.746721,28.736427,28.75036,28.746916,28.831554,28.92691,29.035275,29.11779,29.141161,29.071154,29.068644,29.06333,29.018185,29.034882,28.801441,28.709473,28.692,28.686317,28.76117,28.73792,28.841848,28.761414,28.269253,28.119513,28.170214,28.116533,28.08706,28.177929,28.159948,28.085903,28.100193,28.066307,28.08877,28.138935,28.186434,28.093687,28.005617,28.052176,28.05244,28.131304,28.162014,28.132109,28.109343,28.107265,28.215664,28.246706,28.214806,28.226667,28.19277,28.16817,28.154505,28.053068,28.182304,28.141066,28.099043,28.007734,27.95708,27.922995,27.815264,27.561115,27.467306,27.460306,27.602749,27.676468,27.712864,27.719995,27.74733,27.724382,27.706154,27.743738,27.726295,27.73362,27.799936,27.836857,27.852259,27.845686,27.78033,27.792559,27.84931,27.90373,27.966965,28.01597,28.057632,28.077444,28.089132,28.123755,28.094091,27.977625,27.932549,27.92721,27.865278,27.849676,27.85858,27.856283,27.826319,27.793756,27.78539,27.77959,27.747477,27.685081,27.645636,27.616392,27.589975,27.580963,27.557224,27.52647,27.50037,27.474087,27.46014,27.43747,27.412395,27.452517,27.37899,27.18597,26.838184,26.720562,26.484383,26.352541,26.328293,26.410614,26.396917,26.406412,26.42342,26.42802,26.438559,26.442833,26.442139,26.463963,26.503925,26.46759,26.47291,26.458643,26.45085,26.492603,26.562891,26.499111,26.505167,26.631659,26.866835,26.975391,27.037867,26.96441,26.985113,26.834547,26.281366,25.694977,26.02971,26.087143,25.935766,25.557003,25.663511,24.812897,24.033998,24.53606,24.857018,24.79873,24.154644,24.888182,24.371538,23.361893,23.85227,25.18782,25.636766,25.590132,25.526407,25.365568,25.17262,25.11991,24.984644,24.720207,24.903454,24.927702,24.994286,23.843765,24.366869,24.267948,25.049267,24.9697,23.802027,22.580824,23.875307,24.381643,24.511057,24.282272,23.991169,23.59185,23.316498,23.219011,23.314371,23.524485,24.084549,24.711636,25.151632,25.199446,25.071728,24.856491,24.789858,24.776106,24.724327,24.66991,24.660202,24.61077,24.466454,24.616716,24.67255,24.510115,24.451591,24.602182,24.645512,24.572838,24.391745,24.418058,24.803432,24.8205,24.432495,24.310682,24.399086,24.481874,24.655891,24.822527,24.94012,24.92477,24.950064,24.964104,24.936775,24.857931,24.620853,24.432228,23.7882,23.196724,22.980276,23.499344,23.7011,23.944086,23.81729,23.420933,23.216024,23.020725,22.593302,23.730278,23.735262,23.606092,23.831108,23.687012,23.744865,23.680769,23.583122,23.425955,22.885292,23.079006,22.742191,23.121027,23.267002,23.380539,23.044638,23.178244,22.941412,22.954027,22.747976,22.117897,19.34609,22.02511,21.821556,21.632498,21.962553,21.712317,21.672924,21.82336,22.107584,22.255478,22.090643,21.836092,21.702051,21.853296,21.836723,21.726883,21.588135,21.333323,21.251211,21.226763,21.104645,21.084919,21.186705,21.53116,21.624533,21.409536,21.321165,21.0846,21.063951,21.136482,21.119316,20.966848,20.83051,20.714052,20.507484,20.179502,19.946096,19.997734,20.089523,20.165693,20.238441,20.404331,20.45893,20.41612,20.401049,20.40442,20.418442,20.60399,20.449549,20.211681,20.043262,20.103361,20.888788,20.771658,21.043901,21.69869,21.849403,22.016785,22.093876,22.14757,22.264238,22.34486,22.429605,22.468964,22.476177,22.457312,22.43675,22.419296,22.403725,22.486515,22.654865,22.695894,22.671574,22.470665,22.288929,22.50176,22.69238,22.811527,22.833408,22.830538,22.800293,22.790836,22.769253,22.65044,22.674326,22.5674,22.507198,22.477106,22.515509,22.499195,22.53134,22.5374,22.50533,22.471737,22.443733,22.429012,22.423883,22.414532,22.393265,22.419285,22.441818,22.419264,22.41499,22.425322,22.425175,22.443329,22.4436,22.450813,22.421558,22.396433,22.373611,22.349552,22.357906,22.347273,22.341726,22.346767,22.348608,22.354403,22.36556,22.367722,22.363718,22.390862,22.76037,23.233767,23.348146,23.31802,23.128557,22.83363,22.863016,22.844177,22.823986,22.838854,22.79439,22.767855,22.790861,22.83363,22.820307,22.801924,22.784517,22.73462,22.718443,22.766306,22.717564,22.268011,22.294365,22.149376,22.177471,22.153475, ];

return;
