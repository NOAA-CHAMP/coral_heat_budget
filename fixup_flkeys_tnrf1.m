function stn = fixup_flkeys_tnrf1

  % E.g., http://tds.hycom.org/thredds/dodsC/flkeys.ascii?temperature[0:1463][0][216][258]

  stn = get_flkeys_hycom('tnrf1');

  stn.flkeys_hycom_seatemp.date = datenum(2008,1,1):(6/24):(datenum(2009,1,1)-(0.5/24));
  stn.flkeys_hycom_seatemp.data = [ 23.253332,23.238699,23.072325,23.331615,23.293255,23.522808,23.343931,23.369555,23.22511,22.99628,22.80258,22.245794,22.28128,22.501583,22.254335,22.408628,22.080023,22.63971,22.972906,23.229225,23.329588,23.583502,23.624866,23.781094,23.834646,23.935698,23.989193,24.09083,24.018108,23.867132,23.883066,23.839548,23.721785,23.658297,23.667135,23.650898,23.634558,23.608791,23.580479,23.588673,23.563343,23.566944,23.549097,23.537945,23.551392,23.516302,23.50308,23.534122,23.554264,23.581501,23.57222,23.534023,23.534977,23.51885,23.471525,23.323784,23.2451,23.207176,23.177536,23.155224,23.105368,23.113394,23.120335,23.17576,23.245312,23.369287,23.441917,23.601847,23.675978,23.598885,23.80035,23.865404,23.927721,23.87139,23.86608,23.846977,23.827923,23.814043,23.715763,23.576214,23.428125,23.339205,23.226093,23.235762,23.194933,23.18241,23.18059,23.165602,23.14902,23.130915,23.137503,23.08612,23.097685,23.094742,23.126171,23.277958,23.179794,23.016998,22.86733,22.812065,22.772604,22.75951,22.711555,22.713253,22.769653,22.775692,22.779589,22.699112,22.66739,22.619604,22.561949,22.36915,22.434074,22.280787,22.394005,22.421854,22.392868,22.369238,22.381186,22.363613,22.354612,22.311884,22.317322,22.294407,22.279816,22.24878,22.296375,22.341404,22.411352,22.41124,22.410437,22.381535,23.933973,23.708527,23.71145,23.614601,23.355505,23.247034,23.286875,23.196138,23.005966,22.950123,22.919586,23.000021,22.92247,23.067032,23.278797,23.381775,23.430819,23.649113,23.976156,24.178217,24.536476,24.555643,24.486338,24.560259,24.470621,24.545298,24.578815,24.482662,24.466133,24.370527,24.018356,24.069155,23.850458,23.39938,23.299704,23.270916,23.315033,23.253567,23.199326,23.207436,23.242851,23.260231,23.299198,23.276398,23.251074,23.108715,23.270985,23.993755,24.298481,24.104206,23.935799,23.837534,23.691944,23.456713,23.056612,22.87525,22.937752,22.904512,22.888466,22.901672,22.963055,22.97267,22.848955,22.902416,22.93077,22.85615,22.731838,22.79123,22.680777,22.695526,22.624163,22.569992,22.563835,22.506914,22.53802,22.631548,22.656696,22.59636,22.524925,22.550064,22.717058,22.821888,22.890953,23.745485,24.366348,24.608873,25.004118,25.223875,25.362476,25.606869,25.630657,25.510717,25.280193,25.300684,25.915064,25.529642,25.615692,25.745125,25.521193,25.289795,24.005133,23.549963,23.10178,22.896744,22.315216,22.136559,22.047956,22.025986,21.965055,21.99124,22.056578,22.08745,22.100227,22.104303,22.10276,22.07982,22.061884,22.034998,22.037434,22.086884,22.241926,22.34824,22.582695,22.643642,22.676653,22.785084,22.923101,23.03331,23.028263,23.11652,23.154295,23.09543,23.136839,23.05856,22.717537,22.75304,22.707218,22.678904,22.3024,22.05488,21.883413,21.91192,22.09821,22.28098,22.237656,22.13321,22.073772,22.136705,22.023134,22.019707,21.993053,21.98364,22.006523,22.101963,22.223934,22.381084,22.242805,22.09347,22.009836,22.078987,22.148987,22.173653,22.214031,22.213924,22.329739,22.444904,22.52831,22.732595,22.728285,23.082144,23.402401,23.90994,23.463701,23.823885,23.67549,22.988552,22.816544,22.677397,22.560633,22.611326,22.61967,22.746973,22.795948,22.821758,22.899014,22.882513,22.879822,22.95717,23.063543,23.000954,22.94078,22.907873,22.921354,22.951431,22.999002,23.065666,23.058403,22.989658,23.02935,23.033318,23.051107,23.075266,23.089945,23.28615,23.094316,22.305449,22.09732,22.216806,22.023775,22.304426,22.313412,22.453394,22.468485,22.428623,22.533976,22.649097,22.606672,22.569519,22.587599,22.612898,22.647745,22.678568,22.681885,22.773096,22.85455,22.809303,22.908434,22.885496,22.879143,22.890108,22.977808,23.09412,23.086946,23.134924,23.138819,23.138016,23.153507,23.229704,23.238571,23.222218,23.243975,23.178795,23.124083,23.144018,23.160948,23.152777,23.179964,23.222076,23.344719,23.581024,23.942314,24.049303,24.001125,24.076288,24.118006,24.17097,24.248549,24.32345,24.439365,24.436798,24.341604,24.136786,24.200787,24.276627,24.155565,24.185747,24.23497,24.296616,24.324924,24.337921,24.270052,24.281656,24.317623,24.335094,24.304956,24.323666,24.353886,24.334702,24.39649,24.42491,24.462637,24.535936,24.626762,24.842722,24.703712,24.499128,24.439928,24.506971,24.522566,24.461218,24.371778,24.295658,24.226568,24.184734,23.861494,23.873995,23.870962,23.738682,23.617893,23.702513,23.601564,23.67611,23.71014,23.734707,23.837503,23.907684,23.90119,23.928226,23.932026,23.91569,23.939188,24.094671,24.193605,24.351324,24.353191,24.378094,24.387619,24.185482,24.148518,24.251318,24.52007,24.249813,24.624477,24.57331,24.63037,24.02949,23.925581,24.005377,24.183847,24.49159,24.53061,24.622583,24.820637,24.87122,24.982029,24.949411,24.902592,24.820543,24.798931,24.795418,24.802404,24.839317,24.803738,24.810364,24.840288,24.928247,25.08651,25.089914,25.131914,25.10148,25.067148,25.011047,24.920082,24.85946,24.771988,24.772644,24.758928,24.726137,24.753027,24.751083,24.74561,24.74741,24.756683,24.776707,24.894707,24.98204,25.04565,25.052023,25.054583,25.103796,25.27407,25.377274,25.584648,25.538757,25.1385,25.188763,25.317524,25.317251,25.390059,25.765295,25.409418,25.368801,25.375425,25.326073,25.351772,25.373865,25.368519,25.382334,25.465317,25.774282,26.007668,26.201536,26.309015,26.559616,27.122988,26.873644,27.730515,27.425375,27.160915,26.793116,26.511822,27.665318,26.80018,26.351885,26.8842,26.30546,26.171581,26.043068,26.035837,25.948168,25.884956,25.877241,26.019112,25.986523,26.023027,26.07925,26.207571,26.218832,26.272575,26.377796,26.488976,26.602695,26.9358,26.876904,27.188797,27.652618,27.644451,28.17911,28.232866,28.287762,28.130474,28.188232,27.959496,28.311024,28.176823,28.055561,28.081598,27.94548,27.912462,27.911482,27.839455,27.989103,27.286503,27.419575,27.598795,27.537373,27.395235,27.534946,28.081816,28.535488,28.694576,28.551779,27.928295,27.564009,27.393301,27.478191,27.393534,27.343395,27.268776,27.215408,27.154526,27.14963,27.18245,27.335392,27.442717,27.572426,27.453804,27.371218,27.331488,27.321121,27.287395,27.28117,27.262844,27.196177,27.152817,27.154827,27.154627,27.134113,27.058294,27.04399,27.08539,27.004854,26.86822,26.626421,26.701033,26.78343,26.727102,26.913708,27.139538,27.152088,27.099398,27.098324,27.078255,27.24768,27.22432,27.179493,27.170818,27.181114,27.177086,27.249496,27.259008,27.326557,27.412142,27.47718,27.4986,27.597736,27.629873,27.57492,27.579657,27.586567,27.571547,27.575336,27.56533,27.572674,27.591196,27.606277,27.58103,27.634075,27.714804,27.770689,27.754227,27.729895,27.715206,27.713076,27.699179,27.726906,27.766972,27.823267,27.88923,28.037739,28.063757,28.050943,28.021082,27.979963,27.984253,27.993061,27.956003,27.954786,27.927837,27.944548,27.927727,27.939507,27.919485,27.928629,28.022245,28.117626,28.230825,28.150785,28.067007,28.103811,27.883429,27.87097,28.080404,28.258608,28.265596,28.167545,28.336775,28.482252,28.567469,28.627277,28.59242,28.773964,28.726011,28.92927,29.182364,29.109625,28.931108,29.097952,29.326845,29.24291,28.78091,28.841604,28.90378,28.79896,28.73827,28.713573,28.788683,28.83821,28.7852,28.738169,28.785265,28.746367,28.326988,28.001535,28.051434,28.075874,28.06617,28.198704,28.160545,28.38361,28.482723,28.540869,28.562195,28.568676,28.623402,28.554705,28.504118,28.425312,28.347244,28.267796,28.178417,28.135326,28.079182,27.982658,28.101809,28.264868,28.346527,28.349262,28.332048,28.336048,28.4235,28.479242,28.47769,28.495598,28.504509,28.550207,28.570032,28.542109,28.48655,28.37988,28.312498,28.321222,28.359175,28.275728,28.305595,28.40707,28.425497,28.352123,28.393177,28.421122,28.428894,28.384037,28.356335,28.347567,28.342087,28.35409,28.337612,28.303709,28.357298,28.369907,28.384996,28.377558,28.362217,28.417665,28.387836,28.442244,28.418715,28.385405,28.420866,28.376102,28.214672,28.14433,28.118423,28.0337,28.112373,27.944635,27.99837,28.107332,28.226486,28.168678,28.266632,28.461212,28.51066,28.731796,28.897718,28.79386,28.824842,28.913168,28.913221,28.910719,28.768253,28.73306,28.766726,28.779812,28.631706,28.601748,28.549526,28.599287,28.629827,28.674807,28.628996,28.69106,28.662441,28.519705,28.664186,28.703505,28.824656,28.817595,28.900473,28.91546,28.921825,28.935474,28.973213,29.041899,29.019896,28.997822,28.997955,28.920538,28.898157,28.95709,28.968204,28.964016,29.048853,29.063276,28.957933,28.98747,28.853981,28.750685,28.884705,28.658245,28.592518,28.699852,28.901318,28.900848,28.752563,28.80804,28.892096,29.016031,28.979977,29.126684,29.33403,29.485956,29.46204,29.397831,29.593418,29.657482,29.694344,29.66126,29.64561,29.588755,29.447803,29.315002,29.338757,29.225866,29.176714,29.142076,29.217823,29.238327,29.318817,29.34931,29.384588,29.45772,29.521389,29.514288,29.527748,29.532082,29.526237,29.566338,29.637537,29.670528,29.699043,29.664042,29.620655,29.692488,29.767923,30.093082,29.748735,30.222427,30.395338,30.391174,30.357565,30.517208,30.494907,30.485071,30.460173,30.507984,30.369007,30.35353,30.460375,30.438272,30.455828,30.36569,30.291182,30.329298,30.33438,30.329008,30.559378,30.81969,30.655542,30.688553,30.708881,30.768368,30.710903,30.64396,30.589762,30.498856,30.447453,30.214464,30.10157,30.07471,29.940834,29.807304,29.640867,29.35072,29.171125,29.17527,29.087416,28.849308,28.968225,29.073788,29.003635,29.034506,29.016802,28.998507,28.977966,28.92572,28.891796,28.758833,28.56559,28.824785,28.88736,29.038755,29.19074,29.310675,29.321028,29.39326,29.416353,29.51928,29.51042,29.5002,29.47296,29.44449,29.283865,29.07969,29.225521,29.306864,29.281332,29.369118,29.605711,29.620058,29.742237,29.795729,29.754917,29.699947,29.665575,29.628613,29.560007,29.506788,29.429646,29.356441,29.284803,29.192606,29.097565,29.071299,28.9974,29.015657,29.089458,29.27749,29.382504,29.310535,29.278584,29.262653,29.239151,29.218922,29.2204,29.237734,29.21536,29.204124,29.215754,29.201178,29.214355,29.188286,29.243095,29.23133,29.24035,29.270992,29.28621,29.231064,29.2748,29.274057,29.284615,29.255276,29.18858,29.114975,29.019554,28.840849,28.778006,28.55836,28.371876,28.3606,28.265774,28.273413,28.16833,28.119553,28.130522,28.127457,28.130033,28.156143,28.15147,28.121618,28.068542,28.079103,28.108791,28.096468,28.112173,28.149506,28.168335,28.131649,28.085896,28.0675,28.013456,28.036657,28.072027,28.068626,28.160841,28.320152,28.520159,28.716215,28.841358,28.871758,28.929106,29.20819,29.23283,29.03995,29.01251,28.89552,28.96599,28.881056,28.773266,28.65811,28.602074,28.57403,28.704527,28.730083,28.746502,28.735153,28.70063,28.583193,28.513927,28.371073,28.213205,28.209385,28.14163,28.09764,28.093328,28.11555,28.137302,28.161419,28.145906,28.242073,28.249805,28.146938,28.080788,28.166666,28.213734,28.227654,28.218365,28.190226,28.212744,28.15025,28.151047,28.143232,28.153606,28.11318,28.088219,28.107098,28.081093,28.030561,28.071728,28.09245,28.088924,28.079268,28.054247,28.028107,27.932735,27.918123,27.925985,27.852467,27.85947,27.815348,27.738398,27.79106,27.942993,28.069042,28.101006,28.125465,28.044376,27.984404,27.912325,28.043974,28.018778,27.8749,27.906025,27.927437,27.937101,27.94839,27.976557,27.982706,28.017113,27.964848,27.979147,28.017973,28.009865,27.9984,27.98039,27.934212,27.900475,27.89177,27.91279,27.891502,27.885555,27.886587,27.860527,27.87555,27.902422,27.920412,27.956223,27.958601,27.967257,27.97555,27.911236,27.872124,27.845469,27.81789,27.820168,27.819193,27.828234,27.817354,27.798893,27.77379,27.756264,27.720207,27.703028,27.66456,27.630367,27.61852,27.575201,27.568563,27.529999,27.553307,27.543715,27.58319,27.508913,27.475061,27.500025,27.598063,27.592766,27.566734,27.366917,27.307924,27.293566,27.215698,27.20557,27.193188,27.100391,27.008095,26.989063,26.993423,27.024435,27.015316,26.941557,26.923866,26.926165,26.890835,26.865814,26.867912,26.855173,26.835249,26.838627,26.823322,26.764336,26.770065,26.72021,26.710405,26.730568,26.736616,26.722925,26.722239,26.737179,26.658014,26.61826,26.614712,26.55044,26.43116,26.28813,26.120026,25.876278,25.886726,25.880617,25.86878,25.830452,25.673075,25.793833,25.621292,25.68986,25.70469,25.668688,25.679922,25.606276,25.614813,25.559732,25.657934,25.635834,25.569921,25.292334,25.470335,25.477648,25.543055,25.587484,25.57755,25.597952,25.554234,25.528336,25.521175,25.496876,25.420351,25.330103,25.284237,25.22342,25.368443,25.329859,25.33277,25.276524,25.165413,25.120813,25.101686,25.04236,25.040058,25.010563,24.945482,24.859632,24.85651,24.81338,24.777412,24.792492,24.798721,24.784052,24.743294,24.867413,24.937899,24.940111,24.939365,24.924381,24.921726,24.953989,25.068865,25.082476,25.034084,24.886059,24.898321,24.991096,25.117441,25.123856,25.058737,25.514473,25.97307,25.414234,25.277803,25.38999,25.297768,24.99509,24.923399,24.942585,25.059975,25.000296,24.9299,24.876993,24.776737,24.592884,24.412508,24.424883,24.248371,24.2598,24.118006,24.081724,24.007841,24.03344,23.855581,24.14878,24.131332,24.093275,23.896612,23.773354,23.871479,23.456568,23.60033,23.595863,23.519304,23.576113,23.491861,23.503874,23.433018,23.445543,23.341024,23.416815,23.531996,23.524931,23.556389,23.540354,23.402674,23.156065,23.341248,23.5145,23.244272,23.105268,22.128874,20.404375,20.567446,20.276245,20.143312,20.275274,20.94609,21.104523,21.297249,20.808193,21.09126,21.14757,21.4712,21.351566,21.28617,20.976133,20.9476,20.88511,21.045244,20.981632,20.741575,20.79357,20.951132,20.987638,20.966156,21.094854,21.168549,21.206287,21.437792,21.463945,21.603947,21.869347,22.072975,22.278051,22.592737,22.63403,22.46157,21.672108,21.56485,21.667637,21.168608,21.049377,21.173992,21.309761,21.508827,21.741117,21.79013,21.822014,21.906813,21.93734,22.066809,22.079727,21.972635,22.242859,22.39104,22.267376,22.12906,22.026281,22.326435,22.769176,22.884476,22.7852,22.455809,22.427675,22.379053,22.295223,22.310556,22.33201,22.336721,22.580196,22.751705,22.85279,22.911867,22.686802,22.632532,22.60087,22.551775,22.5279,22.509014,22.512266,22.525545,22.577236,22.746237,22.812876,22.765965,22.725737,22.648968,22.804625,23.153923,23.553778,23.614618,23.54218,23.609756,23.598476,23.598244,23.530895,23.48042,23.426907,23.316704,23.161848,23.061844,22.96562,22.864733,22.856918,22.78758,22.83894,22.810759,22.769943,22.74858,22.723427,22.695335,22.680315,22.661898,22.618559,22.631403,22.652655,22.643625,22.674984,22.669374,22.737328,22.813599,22.811878,22.832582,22.901691,23.199677,23.351273,23.477386,23.486912,23.420147,23.424229,23.264359,23.10034,23.183918,23.176819,23.153084,23.162584,23.204857, ];

return;
