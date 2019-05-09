function stn = fixup_flkeys_plsf1

  % E.g., http://tds.hycom.org/thredds/dodsC/flkeys.ascii?temperature[0:1463][0][209][59]

  stn = get_flkeys_hycom('plsf1');

  stn.flkeys_hycom_seatemp.date = datenum(2008,1,1):(6/24):(datenum(2009,1,1)-(0.5/24));
  stn.flkeys_hycom_seatemp.data = [ 23.671078,23.692911,23.712767,23.711435,23.665041,23.615025,23.53433,23.401644,23.216734,23.022795,22.762627,22.539036,22.39805,22.27333,22.08059,21.802618,21.507345,21.233433,20.829021,20.223555,19.848806,21.673422,22.023975,22.29647,22.552505,22.762362,22.857422,22.914543,22.98672,23.057653,23.029858,22.999687,22.95027,22.856295,22.751913,22.778177,22.820124,22.708431,22.345194,21.927855,21.683807,21.48866,21.47677,21.512247,21.448395,21.46832,21.544697,21.957443,22.30845,22.59445,22.619032,22.484703,22.284128,21.881319,21.79031,21.548326,21.327543,20.932756,21.01409,20.99259,21.08658,21.072105,21.008242,20.980488,21.015385,21.09823,21.188274,21.269793,21.36469,21.46239,21.551659,21.637745,21.690567,21.748562,21.820713,21.867212,21.956419,22.101381,22.261469,22.516903,22.40498,22.28345,22.098925,21.910065,21.766636,21.534275,21.262789,21.085403,21.025896,21.025908,21.114244,21.291794,21.440031,21.563242,21.590033,21.621677,21.648186,21.682436,21.652102,21.595137,21.559343,21.5257,21.537363,21.5545,21.599218,21.660322,21.770744,21.95103,21.805418,21.818192,21.700954,21.595383,21.531769,21.62883,21.75277,21.71883,21.638338,21.565485,21.614101,21.624224,21.546648,21.42223,21.344383,21.25989,21.171717,21.115622,21.091759,21.105904,21.36837,21.516472,21.374367,21.35331,21.396309,21.421936,21.432951,21.461473,21.480677,21.507011,21.485912,21.434183,21.386765,21.34855,21.280209,21.149084,21.031345,21.005798,21.070047,21.086485,21.349556,21.555634,21.449856,21.590063,21.667864,21.873667,21.930923,21.9742,22.050463,22.161568,22.10524,22.060078,21.893127,21.590998,21.615545,21.741507,21.880024,22.023024,22.08303,22.112467,22.118008,22.061794,21.97592,21.86394,21.838856,21.862843,21.885294,21.964876,21.931505,21.91397,21.868177,21.79909,21.769884,21.82096,21.904041,21.928413,21.942226,21.920979,21.874329,21.82974,21.821625,21.824282,21.820583,21.82453,21.837915,21.873558,21.982916,22.10785,22.363762,22.59455,22.567362,22.06374,22.065605,22.062725,22.06503,22.074774,22.086962,22.068748,22.075888,22.095015,22.116861,22.138186,22.174469,22.377165,22.743568,22.829334,23.046026,23.283485,23.171621,23.509558,23.073404,22.887716,23.13702,23.3416,23.60027,23.237581,22.985018,23.198465,23.351915,23.126448,23.243904,23.248194,22.754185,22.259893,22.125954,22.016134,21.884855,21.77695,21.665434,21.55803,21.459543,21.406752,21.388718,21.400589,21.43045,21.487389,21.533215,21.609037,21.638123,21.620523,21.61487,21.591412,21.565632,21.569773,21.555098,21.545038,21.514235,21.506117,21.551142,21.589935,21.71123,21.740122,21.85237,21.987507,21.982916,21.772226,21.72939,21.753649,21.72859,21.836224,21.875534,21.839787,21.846766,21.800785,21.809227,21.792479,21.765953,21.753107,21.750408,21.761024,21.767807,21.781427,21.801832,21.82603,21.841476,21.855892,21.873512,21.913094,21.968138,21.978647,21.994612,21.95686,21.941715,21.93774,22.008692,21.970488,21.958916,21.979824,22.035416,22.185972,22.346638,22.482782,22.672113,23.176151,23.675589,23.89042,23.58148,23.505476,23.025978,22.679428,22.517317,22.452509,22.339407,22.246424,22.29516,22.420784,23.132423,23.392471,23.361677,23.41613,23.451996,23.483322,23.40675,23.028528,22.844135,22.79708,22.772633,22.776058,22.818565,22.831884,22.842804,22.829666,22.88878,23.033167,23.09663,22.917662,22.76126,22.657866,22.524086,22.39782,22.270283,22.204708,22.174416,22.164867,22.181232,22.233294,22.26653,22.308489,22.395906,22.587093,22.667347,22.639511,22.742182,22.87463,22.927675,22.946943,22.996454,23.056778,23.126995,23.163553,23.218533,23.20809,23.211948,23.11342,23.054104,23.085588,23.094753,23.139969,23.178934,23.182657,23.187426,23.198883,23.227215,23.263908,23.300047,23.310823,23.33549,23.353437,23.367552,23.347242,23.316854,23.402184,23.570463,23.687042,23.748766,23.823223,23.92475,23.932852,23.87841,23.884317,23.91002,23.918169,23.902328,24.091963,24.13224,24.07277,24.086117,24.12319,24.152008,24.099695,24.066093,24.102978,24.16232,24.184404,24.247366,24.26814,24.290667,24.284235,24.24497,24.317474,24.364885,24.426064,24.459679,24.538729,24.758728,25.111156,25.04683,25.023746,24.730457,24.643906,24.552065,24.542273,24.450846,24.345083,24.190197,24.02454,23.930044,23.814386,23.674366,23.597912,23.625006,23.671396,23.689007,23.734365,23.77478,23.80198,23.816109,23.82434,23.812046,23.813515,23.848656,24.020382,24.074396,24.039425,24.080202,24.221216,24.433546,24.358633,24.468592,24.598896,24.447685,24.435392,24.30723,24.186052,24.245567,24.179777,23.939646,23.789387,23.765005,23.799328,23.806633,23.85347,23.896708,23.857174,23.805126,23.827847,23.772587,23.729622,23.66997,23.660923,23.742544,23.753962,23.653069,23.57336,23.57812,23.543009,23.5547,23.622648,23.756313,23.792383,23.83366,23.869913,23.900011,23.892746,23.89107,23.854923,23.804937,23.801445,23.688297,23.642769,23.543673,23.365444,23.501865,23.725912,23.904661,24.09718,24.261372,24.301931,24.289429,24.31752,24.39745,24.509289,24.595259,24.667961,24.80605,25.001244,25.056202,24.858261,24.67172,24.726017,24.816708,24.878649,24.967722,25.014582,24.972483,25.06318,25.0615,25.046434,25.187332,25.695017,25.943361,26.101284,26.333435,26.609192,26.636427,26.652905,26.846664,27.135267,27.07621,27.025167,27.063858,27.084528,27.198484,26.973253,26.615065,26.70405,26.562487,26.302696,26.348541,26.050125,25.638348,25.361523,25.116999,25.048191,24.973946,24.804525,24.777678,24.925224,24.934027,25.038368,25.28559,25.526384,25.86008,26.071754,25.989595,26.272215,26.46827,26.325493,26.378538,26.415854,26.265535,26.360733,26.453094,26.543827,26.530153,26.230078,26.1178,26.268847,26.39202,26.39402,26.468157,26.635569,26.703886,26.774933,26.83945,26.824684,26.695433,26.805624,27.052483,27.093313,26.951988,27.214506,27.317488,27.247906,27.221777,27.308872,27.104786,26.906996,26.630627,26.270763,26.020287,25.681768,24.869576,24.675726,25.867996,26.272074,26.413845,26.42644,26.425085,26.361649,26.34166,26.292046,26.31657,26.388268,26.45532,26.485939,26.538683,26.539413,26.594706,26.674547,26.899223,27.059063,27.087755,27.248383,27.417053,27.422573,27.394997,27.400606,27.350267,27.39942,27.564558,27.606318,27.594505,27.357033,27.326141,27.403334,27.471775,27.109709,27.01048,26.92042,26.98413,26.84401,26.92846,26.777824,26.46632,25.723679,25.582504,25.831734,26.135904,26.36499,26.478191,26.55909,26.698853,26.531668,26.488873,26.672174,26.62279,26.69038,27.097748,27.327906,27.37708,27.407633,27.494078,27.630135,27.648832,27.560337,27.575537,27.596571,27.633799,27.619799,27.57961,27.559757,27.583944,27.599218,27.58384,27.583853,27.530317,27.534687,27.581247,27.629404,27.646273,27.655506,27.674778,27.71005,27.812492,27.8551,27.91509,27.95779,28.049597,28.126005,28.051048,27.978477,27.953457,27.993948,27.99565,27.945438,27.935982,28.081474,28.13572,28.122044,28.161304,28.251604,28.233156,28.302647,28.338898,28.272999,28.495623,28.769753,28.970201,29.02034,28.8549,28.914309,28.94818,28.640766,28.483246,28.476933,28.50095,28.51565,28.502056,28.526363,28.512636,28.482824,28.343668,28.221876,27.997204,27.892097,27.82061,27.822504,27.881691,27.842184,27.89145,27.94386,27.988777,28.036818,28.1321,28.238916,28.328989,28.3839,28.405602,28.537924,28.654737,28.655146,28.620602,28.541338,28.516142,28.480467,28.488783,28.459044,28.403927,28.352991,28.30801,28.349453,28.375334,28.415838,28.404415,28.431732,28.463686,28.503382,28.553574,28.580532,28.570068,28.578857,28.575087,28.54963,28.537056,28.501188,28.451618,28.452137,28.470003,28.472904,28.477266,28.46426,28.450617,28.465403,28.474695,28.539577,28.519117,28.507967,28.52903,28.548334,28.519514,28.454506,28.408342,28.423628,28.469643,28.506708,28.537516,28.564207,28.595043,28.618607,28.632969,28.6427,28.669876,28.736256,28.817516,28.861446,28.784899,28.721188,28.690008,28.662632,28.6066,28.595293,28.630949,28.623205,28.639675,28.70346,28.842459,29.046576,29.267494,29.316519,29.096468,29.025532,29.03484,28.972895,28.95144,28.936146,28.933104,28.891066,28.880749,28.85353,28.872574,28.849882,28.799072,28.792109,28.815142,28.856514,28.861929,28.878286,28.908546,28.934998,28.917936,28.905338,28.923458,28.93949,28.953894,28.970215,28.990963,28.997255,28.989944,29.015617,29.053553,29.003681,29.000319,29.024754,29.100222,29.261488,29.241978,29.328285,29.404459,29.465235,29.488302,29.470766,29.416235,29.395084,29.368048,29.31717,29.365473,29.361448,29.361607,29.375208,29.383917,29.409105,29.365597,29.411964,29.48523,29.51401,29.480175,29.481153,29.519888,29.64405,29.681189,29.770134,29.853374,29.789799,29.717054,29.659805,29.734255,29.665834,29.623629,29.585253,29.613445,29.61422,29.554962,29.62296,29.647251,29.621979,29.601353,29.595636,29.623962,29.638363,29.651217,29.760902,29.84845,29.803444,29.723341,29.691658,29.731297,29.697626,29.760836,29.764812,29.77046,29.853933,29.917868,29.91311,29.894375,29.907982,29.901089,29.929152,29.933,29.969828,29.945518,29.947563,29.939411,30.028145,30.118511,30.111876,30.08488,30.145847,30.11095,30.096985,30.142088,30.347408,30.354696,30.283804,30.369574,30.396818,30.302954,30.278893,30.27063,30.1976,30.164196,30.127375,30.025457,29.952858,29.851194,29.59996,29.341291,29.206564,29.26922,29.4308,29.388186,29.427204,29.42286,29.412766,29.41772,29.385475,29.36862,29.387476,29.3638,29.32063,29.38334,29.353273,29.319077,29.31501,29.361645,29.436321,29.531727,29.53811,29.48286,29.47529,29.562916,29.582716,29.614023,29.607056,29.603252,29.580486,29.596071,29.60018,29.64625,29.67567,29.694206,29.765045,29.685505,29.635712,29.488523,29.371206,29.41118,29.540022,29.589762,29.59522,29.517292,29.454462,29.09242,28.758333,28.418318,28.559923,28.800488,28.667194,28.617058,28.640444,28.719828,28.778831,28.987373,28.930996,28.921608,28.98785,29.00389,28.978743,28.933752,28.881748,28.815153,28.820597,28.805387,28.812244,28.840576,28.834745,28.860617,28.876339,28.849123,28.814009,28.819868,28.66099,28.40429,28.33939,28.332367,28.319424,28.336435,28.449492,28.443811,28.392448,28.476582,28.461588,28.620827,29.483429,29.173962,28.27563,28.152615,28.058144,28.110546,28.136147,28.183773,28.292395,28.486706,28.80747,28.867949,28.870016,28.792295,28.745518,28.702948,28.685194,28.726147,28.776873,28.840021,28.77608,28.782993,28.764751,28.810038,28.860847,28.897388,28.912039,28.899324,28.963057,29.03395,29.054186,29.152443,29.27975,29.239752,29.22262,29.216543,29.1865,29.110254,28.983805,28.861952,28.663832,28.557352,28.476244,28.358524,28.245724,28.204657,28.192291,28.207224,28.208395,28.239384,28.270617,28.285814,28.305733,28.41108,28.439009,28.445684,28.401108,28.373066,28.359297,28.24328,28.143702,28.014925,27.989021,27.920753,27.892315,27.968403,28.05587,28.13121,28.222466,28.324827,28.309568,28.31397,28.34969,28.344976,28.23031,28.132195,28.154654,28.08131,27.883795,27.882252,27.944792,27.999565,28.003681,27.97569,27.984518,28.005653,28.055065,28.081223,28.087261,28.100132,28.143082,28.162113,28.164253,28.212177,28.219124,28.232542,28.157938,28.1594,28.12384,28.09788,28.09307,28.074492,28.044895,28.019289,27.996767,27.967789,27.964752,27.965801,27.994226,28.011623,28.024097,28.0359,28.04824,28.073874,28.10542,28.11434,28.120016,28.135956,28.145344,28.131557,28.115078,28.101189,28.082844,28.09153,28.11184,28.167082,28.181553,28.187662,28.164995,28.103369,28.072638,28.026,28.008339,27.97546,27.938145,27.91512,27.91168,27.936354,27.952066,27.955534,27.954643,27.943176,27.931633,27.928497,27.908009,27.843332,27.723833,27.582823,27.374226,27.262085,27.322777,27.332592,27.383457,27.445107,27.412092,27.31503,27.207615,27.148354,27.064457,26.994698,27.01243,27.037205,27.085249,27.162945,27.201488,27.170883,27.128069,27.02497,26.871746,26.715202,26.593287,26.456444,26.316591,26.15452,26.076277,26.251713,26.45095,26.417337,26.768179,26.791634,26.793142,26.09185,26.030542,26.078226,26.147217,26.124329,26.421322,26.273859,25.974077,26.053814,26.33264,26.58448,26.542757,26.493675,26.458956,26.476006,26.442007,26.38559,26.314575,26.240791,26.195435,26.165358,26.1465,26.104666,26.058191,25.987715,25.894302,25.77995,25.687872,25.582863,25.476963,25.369606,25.290054,25.225815,25.159678,25.081978,25.061789,25.057531,25.079584,25.065283,25.020943,25.001461,25.007877,25.03732,25.02986,25.016556,25.017899,25.019714,25.05711,25.074514,25.103262,25.124523,25.133614,25.11889,25.08052,25.029308,24.988441,24.969402,24.97228,24.996881,25.028397,25.065369,25.102695,25.1335,25.16845,25.211039,25.246038,25.238064,25.169855,25.115896,25.063997,24.967623,24.869543,24.793398,24.724642,24.668894,24.609003,24.582338,24.561646,24.539076,24.548574,24.594166,24.658611,24.719082,24.766022,24.802006,24.843815,24.880358,24.939758,24.941265,24.94808,24.929829,24.878004,24.80109,24.732618,24.645435,24.564695,24.528265,24.503735,24.500832,24.50068,24.523935,24.518692,24.44456,24.363272,24.251528,24.13248,24.034702,23.96022,23.93302,23.902487,23.90039,23.912079,23.918545,23.914112,23.86477,23.78215,23.682302,23.569143,23.439617,23.327463,23.223366,23.126682,23.032234,22.962221,22.942366,22.98597,23.068022,23.18629,23.30316,23.42745,23.4855,23.502644,23.48894,23.47545,23.444714,23.41629,23.390089,23.370916,23.349468,23.391008,23.394655,23.323984,23.317287,23.29725,23.339588,23.339436,23.335463,23.338257,23.32921,23.317387,23.27259,23.20761,23.187086,23.25173,23.289999,23.267048,23.13761,22.977476,22.844463,22.73205,22.6258,22.500214,22.33747,22.187588,22.09802,22.034002,21.96503,21.922174,21.89266,21.861286,21.817204,21.839659,21.876324,21.903908,21.88342,21.950018,22.043062,22.11825,22.165154,22.205713,22.244848,22.251255,22.216494,22.19844,22.212605,22.22205,22.251741,22.28995,22.320786,22.331821,22.345966,22.39812,22.489529,22.524302,22.529558,22.53445,22.468908,22.442093,22.367725,22.261456,22.142427,22.030315,21.867956,21.749025,21.732895,21.766832,21.79483,21.945646,22.167141,22.279482,22.29128,22.292639,22.326342,22.376728,22.37199,22.432755,22.38246,22.414253,22.367296,22.36295,22.358892,22.356596,22.3368,22.336006,22.355518,22.37053,22.381874,22.427137,22.432047,22.370346,22.343435,22.312645,22.255117,22.198803,22.14472,22.163046,22.09597,21.782822,21.95072,21.795967,21.337042,21.177765,21.295551,21.383736,21.50098,21.685476,21.837812,21.992167,22.114185,22.160439,22.172178,22.139297,22.080061,22.034561,21.95824,21.860538,21.719908,21.570187,21.50888,21.479822,21.48217,21.480768,21.505507,21.67004,21.780067,21.902325,21.950382,22.067865,22.193867,22.21557,22.200048,22.32827, ];

return;
