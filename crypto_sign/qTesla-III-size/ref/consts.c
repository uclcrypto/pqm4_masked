/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: constants for the NTT
**************************************************************************************/

#include <stdint.h>
#include "params.h"
#include "poly.h"

const poly zeta ={

671800, 4181291, 975654, 970836, 1459996, 2949013, 1578790, 3375131, 177347, 2024971, 3299069, 2879655, 1061156, 3772041, 1726661, 2646527, 
224962, 3106510, 1764167, 3790159, 110295, 277183, 2296602, 1995237, 1574725, 1473236, 1081285, 144829, 114244, 719647, 4114328, 917441, 
4188270, 3805772, 261389, 52393, 2185303, 1021265, 2167874, 2986441, 3886274, 2191966, 284211, 3446813, 1389427, 2107810, 1173125, 1597161, 
3753261, 1373052, 793684, 4091628, 1677907, 4164049, 1948749, 2758369, 1027640, 1118203, 891820, 1309242, 1810791, 1863364, 2587868, 1541007, 
4104068, 675426, 1402433, 2557508, 1068970, 1940808, 3957823, 798456, 4092960, 3262467, 1793460, 658044, 1978921, 1367494, 3136736, 2360480, 
941550, 37800, 1919065, 3032526, 581001, 3323192, 299785, 3114533, 545048, 2845265, 1891473, 102035, 2256179, 221259, 1796623, 504470, 
377401, 3184337, 3107383, 606431, 200460, 3770995, 986925, 207500, 3712747, 1696453, 4158053, 3530443, 32005, 3222743, 3918763, 3574153, 
2768592, 2608835, 1856937, 905294, 214652, 4154226, 2876170, 2651799, 1098009, 3905542, 3763042, 3055325, 1438567, 969841, 2397140, 3637385, 
3779810, 644984, 1638607, 498549, 3404792, 4055115, 9472, 315805, 1796876, 972163, 3025826, 3334639, 2290368, 2552107, 160996, 3282568, 
2279239, 1305163, 2304247, 603598, 3803059, 2582009, 3202587, 1094032, 1195417, 2879417, 1648902, 542294, 3085586, 3325229, 4177450, 150226, 
890698, 503530, 3122945, 1929018, 3309179, 1075767, 2185016, 276011, 1620579, 1349757, 454010, 3835301, 3658519, 2369797, 203221, 2116132, 
1371940, 3499700, 2991078, 3597638, 942280, 506271, 701979, 1853372, 2165162, 2830558, 2083508, 3582128, 4177826, 2623861, 3740436, 725559, 
791017, 595361, 2192451, 878351, 1919935, 1730363, 165115, 3011415, 539166, 4049306, 2512830, 3633034, 3743092, 1721797, 356766, 3860922, 
551806, 520752, 1492250, 3020875, 296084, 3951086, 3702654, 1541222, 2760082, 2967699, 1811892, 1913148, 3121111, 2583448, 37791, 2289197, 
228811, 315449, 2711375, 2035264, 998876, 684125, 1377229, 1723513, 2093137, 1181754, 3978572, 1168111, 1295590, 1870157, 2992279, 1610031, 
2052968, 3195982, 3195020, 2498826, 2430997, 1447298, 2178224, 1573739, 3444420, 2425537, 3066466, 1895376, 3494178, 2341084, 2603206, 2810264, 
3665075, 4030046, 232980, 3770527, 2425457, 3193512, 1906687, 3838549, 1081341, 3385499, 343154, 3648238, 3066045, 3502707, 903006, 2216085, 
2477447, 3769256, 2700907, 2899931, 3094342, 404354, 2325640, 4161594, 1153616, 2601633, 624385, 56418, 4122920, 303574, 3474524, 3047326, 
3446806, 1755473, 1289687, 3030484, 745529, 2037059, 1126174, 3508536, 3263841, 1057863, 2424516, 3666380, 2238799, 1918076, 1096624, 666757, 
2414037, 4105141, 86489, 236751, 2175830, 2842379, 3751432, 366978, 1727916, 627613, 2576775, 231383, 2352896, 1039386, 650148, 3849095, 
2893195, 2813545, 2172937, 1389954, 1261168, 1470030, 832830, 3548304, 2585258, 3650945, 3733752, 797947, 4183412, 261772, 374082, 3717015, 
1306771, 591941, 3320862, 3969254, 3730288, 4153963, 2641916, 706453, 3574687, 687011, 3723863, 518936, 674472, 2242626, 174183, 3560884, 
3969544, 425417, 789235, 4183047, 4027225, 982625, 2075760, 2392513, 2538340, 3022462, 1997528, 356548, 3730142, 1536313, 1202696, 1344848, 
3103217, 2383022, 1762142, 2994989, 3102783, 3072599, 1517632, 2024436, 2534641, 147328, 2356097, 190578, 2587663, 440306, 2374767, 3182600, 
680532, 484370, 4131095, 3009332, 3562207, 976019, 3613316, 3033006, 3743622, 4136404, 1605237, 66645, 3859240, 908865, 4051121, 2726336, 
3637443, 2340134, 813357, 3985220, 2868243, 3650243, 1684957, 3023114, 2402323, 1820096, 1764462, 1049670, 2260628, 4976, 3760346, 3157996, 
3573461, 1006628, 454916, 4159906, 337885, 22277, 520578, 2607705, 2561874, 503606, 1415232, 3823408, 3829828, 554510, 914738, 3838536, 
653901, 664244, 3918457, 361056, 515834, 2583400, 2666144, 1562200, 2635470, 3523620, 2847787, 281762, 1416774, 4047010, 2739024, 1492985, 
2613083, 2116726, 4076288, 4141191, 3357856, 741301, 977038, 4028938, 1661277, 2769449, 3571042, 2601104, 57237, 3026729, 3478919, 87366, 
3697654, 2676961, 3932341, 2883942, 3200147, 623723, 871365, 763732, 2354543, 661482, 1442350, 148821, 966821, 2154509, 1229800, 2252524, 
1712762, 687319, 1231124, 225814, 127675, 2786959, 2996601, 1997279, 1410197, 2759369, 254896, 2633749, 743622, 2420984, 594581, 1359068, 
3724994, 3338166, 473524, 3323698, 110693, 2630130, 3742099, 1392129, 1263087, 1474128, 964094, 3338617, 2682625, 2350723, 4039051, 2437147, 
2003303, 1029372, 3710675, 1198388, 4047402, 337401, 959139, 3673320, 3269977, 1757086, 3846011, 3052386, 1555886, 1213798, 1730449, 574426, 
3730903, 4058825, 3075, 3232877, 597243, 584901, 3208277, 423060, 3216342, 3727213, 89571, 709528, 3722455, 112585, 4199553, 578587, 
1727014, 3010665, 1118724, 3088559, 458307, 695931, 2551953, 3462204, 2654347, 2908501, 3034211, 3511237, 3734268, 2443875, 270514, 776347, 
683036, 1526569, 521044, 2352920, 557737, 4056083, 2391161, 2389563, 2293979, 2581739, 2738173, 2545480, 1008072, 3577574, 1673061, 4116273, 
133058, 1222352, 1144238, 882222, 3000625, 4046931, 141504, 1904001, 1035854, 3807884, 2398461, 446181, 2041489, 1148183, 2291458, 3675915, 
255124, 369448, 3016249, 2025225, 3237403, 2220199, 3134791, 2587255, 3220754, 3366174, 132697, 3383227, 1358468, 1158291, 2321651, 2869559, 
1425523, 4054733, 69091, 3521561, 2453355, 2968118, 2968833, 3185424, 988606, 3025251, 1154802, 24092, 1305476, 3938667, 3405455, 2280837, 
2987149, 1576181, 3812113, 481232, 1911887, 2305037, 3637072, 515558, 3183843, 3460525, 2134536, 3376047, 849276, 912675, 3126131, 3349335, 
1736653, 247313, 307171, 2906949, 2483567, 3951115, 449581, 3211241, 119780, 3050685, 312715, 129516, 1413964, 3626707, 1834389, 739674, 
2166987, 1898439, 3247386, 543470, 3893129, 1952324, 4010533, 2663329, 1611039, 4159354, 3221090, 4011118, 456104, 4128401, 3481956, 1341852, 
1346376, 1373597, 1886912, 2289124, 2035164, 3802432, 4020200, 1440583, 131860, 2447356, 1147783, 3884191, 36600, 1417517, 3115113, 4106357, 
2209232, 3913295, 2079509, 2915453, 253356, 2093028, 3105753, 420898, 3641863, 2237777, 589597, 3471638, 1556385, 1574364, 2961455, 2414774, 
2532838, 3894119, 2561579, 1825751, 2610770, 4095615, 2366084, 1696032, 2935352, 1982899, 3940806, 962691, 2874348, 2295425, 3088987, 1724605, 
138760, 2611152, 2321223, 3862854, 977071, 3373271, 2119442, 2444640, 184156, 2401204, 2250096, 3883423, 24318, 1799015, 2709027, 2477092, 
2937887, 872546, 348220, 49520, 266109, 1166709, 1470353, 712356, 4162049, 2520023, 1093919, 3371334, 1529777, 3549597, 3033168, 3626405, 
317815, 972428, 3325840, 1416192, 1615043, 3225312, 49030, 591050, 3470933, 533400, 905783, 2128579, 2589779, 1556207, 3295501, 3128246, 
1323037, 2836289, 1222103, 2635029, 764092, 1785154, 1271391, 407326, 3293361, 697832, 1957938, 26925, 1909470, 921060, 1189793, 452905, 
177180, 3986522, 3612073, 2634482, 3811697, 2155464, 3184049, 3773906, 4155559, 890604, 965647, 946702, 2980153, 2514794, 1634712, 1413135, 
2059115, 1095948, 1094602, 3286386, 1617289, 2234906, 2942756, 2831603, 2790364, 4110867, 3267277, 2818115, 997189, 3975212, 1919457, 3294782, 
42631, 3979780, 677949, 1074671, 4007873, 2013224, 4064265, 1404667, 1266413, 1753048, 1480954, 2251688, 3671233, 3337348, 3023835, 704482, 
1867102, 2290506, 1202801, 3686892, 3618479, 3297031, 2477670, 3415258, 2889498, 2106378, 361488, 1478812, 3536666, 645275, 2793501, 2983604, 
3150760, 1136423, 2629214, 3144871, 1095947, 2432448, 3144124, 1562104, 3685583, 2519659, 1745378, 275993, 1028739, 4053547, 3139341, 791685, 
205316, 940435, 3044250, 3537550, 2347550, 2748749, 216515, 2376693, 3994272, 758809, 336837, 4138282, 254982, 2087732, 1443586, 2090448, 
2407213, 2192231, 584225, 1528366, 714102, 2781015, 1061159, 144894, 2251444, 1706143, 3064185, 1082774, 1212561, 1964667, 1808852, 1281436, 
380192, 1938622, 3594224, 2865093, 1814198, 3709791, 3557452, 641073, 3449310, 2797672, 1886229, 2374072, 1947652, 910530, 4110612, 3688785, 
2761424, 2192378, 1210992, 432423, 3990493, 3710041, 3364266, 1402625, 1430941, 466915, 2307343, 3969361, 3041855, 1636011, 2336989, 4083954, 
1752367, 1468975, 4003767, 2752277, 2639144, 2435428, 168292, 2731409, 2173963, 313121, 1885409, 1792411, 105750, 1595875, 3535511, 1917121, 
3348968, 3600516, 3025874, 1234611, 387230, 254793, 2267177, 423073, 3643782, 2241875, 720861, 3996710, 2066073, 1031892, 2436415, 3356685, 
3628852, 1131491, 315588, 3085726, 4060906, 3713538, 561022, 142143, 137017, 4091465, 525060, 523088, 2581256, 2546361, 529201, 1724592, 
3917913, 4096490, 1689933, 575672, 2633453, 2453964, 3882580, 236313, 394169, 2731312, 3191196, 135139, 1208112, 2180950, 4051722, 330078, 
4161293, 3314132, 1075088, 3797989, 958522, 1974573, 3610471, 3368492, 629863, 3712506, 281606, 4189621, 1437509, 2515187, 1936773, 3150875, 
797596, 4050969, 2506561, 2023050, 3235484, 2216101, 3003527, 569898, 2081018, 3678953, 3392925, 857476, 1224594, 2996526, 3160227, 35843, 
};
const poly zetainv ={

1046366, 1210067, 2981999, 3349117, 813668, 527640, 2125575, 3636695, 1203066, 1990492, 971109, 2183543, 1700032, 155624, 3408997, 1055718, 
2269820, 1691406, 2769084, 16972, 3924987, 494087, 3576730, 838101, 596122, 2232020, 3248071, 408604, 3131505, 892461, 45300, 3876515, 
154871, 2025643, 2998481, 4071454, 1015397, 1475281, 3812424, 3970280, 324013, 1752629, 1573140, 3630921, 2516660, 110103, 288680, 2482001, 
3677392, 1660232, 1625337, 3683505, 3681533, 115128, 4069576, 4064450, 3645571, 493055, 145687, 1120867, 3891005, 3075102, 577741, 849908, 
1770178, 3174701, 2140520, 209883, 3485732, 1964718, 562811, 3783520, 1939416, 3951800, 3819363, 2971982, 1180719, 606077, 857625, 2289472, 
671082, 2610718, 4100843, 2414182, 2321184, 3893472, 2032630, 1475184, 4038301, 1771165, 1567449, 1454316, 202826, 2737618, 2454226, 122639, 
1869604, 2570582, 1164738, 237232, 1899250, 3739678, 2775652, 2803968, 842327, 496552, 216100, 3774170, 2995601, 2014215, 1445169, 517808, 
95981, 3296063, 2258941, 1832521, 2320364, 1408921, 757283, 3565520, 649141, 496802, 2392395, 1341500, 612369, 2267971, 3826401, 2925157, 
2397741, 2241926, 2994032, 3123819, 1142408, 2500450, 1955149, 4061699, 3145434, 1425578, 3492491, 2678227, 3622368, 2014362, 1799380, 2116145, 
2763007, 2118861, 3951611, 68311, 3869756, 3447784, 212321, 1829900, 3990078, 1457844, 1859043, 669043, 1162343, 3266158, 4001277, 3414908, 
1067252, 153046, 3177854, 3930600, 2461215, 1686934, 521010, 2644489, 1062469, 1774145, 3110646, 1061722, 1577379, 3070170, 1055833, 1222989, 
1413092, 3561318, 669927, 2727781, 3845105, 2100215, 1317095, 791335, 1728923, 909562, 588114, 519701, 3003792, 1916087, 2339491, 3502111, 
1182758, 869245, 535360, 1954905, 2725639, 2453545, 2940180, 2801926, 142328, 2193369, 198720, 3131922, 3528644, 226813, 4163962, 911811, 
2287136, 231381, 3209404, 1388478, 939316, 95726, 1416229, 1374990, 1263837, 1971687, 2589304, 920207, 3111991, 3110645, 2147478, 2793458, 
2571881, 1691799, 1226440, 3259891, 3240946, 3315989, 51034, 432687, 1022544, 2051129, 394896, 1572111, 594520, 220071, 4029413, 3753688, 
3016800, 3285533, 2297123, 4179668, 2248655, 3508761, 913232, 3799267, 2935202, 2421439, 3442501, 1571564, 2984490, 1370304, 2883556, 1078347, 
911092, 2650386, 1616814, 2078014, 3300810, 3673193, 735660, 3615543, 4157563, 981281, 2591550, 2790401, 880753, 3234165, 3888778, 580188, 
1173425, 656996, 2676816, 835259, 3112674, 1686570, 44544, 3494237, 2736240, 3039884, 3940484, 4157073, 3858373, 3334047, 1268706, 1729501, 
1497566, 2407578, 4182275, 323170, 1956497, 1805389, 4022437, 1761953, 2087151, 833322, 3229522, 343739, 1885370, 1595441, 4067833, 2481988, 
1117606, 1911168, 1332245, 3243902, 265787, 2223694, 1271241, 2510561, 1840509, 110978, 1595823, 2380842, 1645014, 312474, 1673755, 1791819, 
1245138, 2632229, 2650208, 734955, 3616996, 1968816, 564730, 3785695, 1100840, 2113565, 3953237, 1291140, 2127084, 293298, 1997361, 100236, 
1091480, 2789076, 4169993, 322402, 3058810, 1759237, 4074733, 2766010, 186393, 404161, 2171429, 1917469, 2319681, 2832996, 2860217, 2864741, 
724637, 78192, 3750489, 195475, 985503, 47239, 2595554, 1543264, 196060, 2254269, 313464, 3663123, 959207, 2308154, 2039606, 3466919, 
2372204, 579886, 2792629, 4077077, 3893878, 1155908, 4086813, 995352, 3757012, 255478, 1723026, 1299644, 3899422, 3959280, 2469940, 857258, 
1080462, 3293918, 3357317, 830546, 2072057, 746068, 1022750, 3691035, 569521, 1901556, 2294706, 3725361, 394480, 2630412, 1219444, 1925756, 
801138, 267926, 2901117, 4182501, 3051791, 1181342, 3217987, 1021169, 1237760, 1238475, 1753238, 685032, 4137502, 151860, 2781070, 1337034, 
1884942, 3048302, 2848125, 823366, 4073896, 840419, 985839, 1619338, 1071802, 1986394, 969190, 2181368, 1190344, 3837145, 3951469, 530678, 
1915135, 3058410, 2165104, 3760412, 1808132, 398709, 3170739, 2302592, 4065089, 159662, 1205968, 3324371, 3062355, 2984241, 4073535, 90320, 
2533532, 629019, 3198521, 1661113, 1468420, 1624854, 1912614, 1817030, 1815432, 150510, 3648856, 1853673, 3685549, 2680024, 3523557, 3430246, 
3936079, 1762718, 472325, 695356, 1172382, 1298092, 1552246, 744389, 1654640, 3510662, 3748286, 1118034, 3087869, 1195928, 2479579, 3628006, 
7040, 4094008, 484138, 3497065, 4117022, 479380, 990251, 3783533, 998316, 3621692, 3609350, 973716, 4203518, 147768, 475690, 3632167, 
2476144, 2992795, 2650707, 1154207, 360582, 2449507, 936616, 533273, 3247454, 3869192, 159191, 3008205, 495918, 3177221, 2203290, 1769446, 
167542, 1855870, 1523968, 867976, 3242499, 2732465, 2943506, 2814464, 464494, 1576463, 4095900, 882895, 3733069, 868427, 481599, 2847525, 
3612012, 1785609, 3462971, 1572844, 3951697, 1447224, 2796396, 2209314, 1209992, 1419634, 4078918, 3980779, 2975469, 3519274, 2493831, 1954069, 
2976793, 2052084, 3239772, 4057772, 2764243, 3545111, 1852050, 3442861, 3335228, 3582870, 1006446, 1322651, 274252, 1529632, 508939, 4119227, 
727674, 1179864, 4149356, 1605489, 635551, 1437144, 2545316, 177655, 3229555, 3465292, 848737, 65402, 130305, 2089867, 1593510, 2713608, 
1467569, 159583, 2789819, 3924831, 1358806, 682973, 1571123, 2644393, 1540449, 1623193, 3690759, 3845537, 288136, 3542349, 3552692, 368057, 
3291855, 3652083, 376765, 383185, 2791361, 3702987, 1644719, 1598888, 3686015, 4184316, 3868708, 46687, 3751677, 3199965, 633132, 1048597, 
446247, 4201617, 1945965, 3156923, 2442131, 2386497, 1804270, 1183479, 2521636, 556350, 1338350, 221373, 3393236, 1866459, 569150, 1480257, 
155472, 3297728, 347353, 4139948, 2601356, 70189, 462971, 1173587, 593277, 3230574, 644386, 1197261, 75498, 3722223, 3526061, 1023993, 
1831826, 3766287, 1618930, 4016015, 1850496, 4059265, 1671952, 2182157, 2688961, 1133994, 1103810, 1211604, 2444451, 1823571, 1103376, 2861745, 
3003897, 2670280, 476451, 3850045, 2209065, 1184131, 1668253, 1814080, 2130833, 3223968, 179368, 23546, 3417358, 3781176, 237049, 645709, 
4032410, 1963967, 3532121, 3687657, 482730, 3519582, 631906, 3500140, 1564677, 52630, 476305, 237339, 885731, 3614652, 2899822, 489578, 
3832511, 3944821, 23181, 3408646, 472841, 555648, 1621335, 658289, 3373763, 2736563, 2945425, 2816639, 2033656, 1393048, 1313398, 357498, 
3556445, 3167207, 1853697, 3975210, 1629818, 3578980, 2478677, 3839615, 455161, 1364214, 2030763, 3969842, 4120104, 101452, 1792556, 3539836, 
3109969, 2288517, 1967794, 540213, 1782077, 3148730, 942752, 698057, 3080419, 2169534, 3461064, 1176109, 2916906, 2451120, 759787, 1159267, 
732069, 3903019, 83673, 4150175, 3582208, 1604960, 3052977, 44999, 1880953, 3802239, 1112251, 1306662, 1505686, 437337, 1729146, 1990508, 
3303587, 703886, 1140548, 558355, 3863439, 821094, 3125252, 368044, 2299906, 1013081, 1781136, 436066, 3973613, 176547, 541518, 1396329, 
1603387, 1865509, 712415, 2311217, 1140127, 1781056, 762173, 2632854, 2028369, 2759295, 1775596, 1707767, 1011573, 1010611, 2153625, 2596562, 
1214314, 2336436, 2911003, 3038482, 228021, 3024839, 2113456, 2483080, 2829364, 3522468, 3207717, 2171329, 1495218, 3891144, 3977782, 1917396, 
4168802, 1623145, 1085482, 2293445, 2394701, 1238894, 1446511, 2665371, 503939, 255507, 3910509, 1185718, 2714343, 3685841, 3654787, 345671, 
3849827, 2484796, 463501, 573559, 1693763, 157287, 3667427, 1195178, 4041478, 2476230, 2286658, 3328242, 2014142, 3611232, 3415576, 3481034, 
466157, 1582732, 28767, 624465, 2123085, 1376035, 2041431, 2353221, 3504614, 3700322, 3264313, 608955, 1215515, 706893, 2834653, 2090461, 
4003372, 1836796, 548074, 371292, 3752583, 2856836, 2586014, 3930582, 2021577, 3130826, 897414, 2277575, 1083648, 3703063, 3315895, 4056367, 
29143, 881364, 1121007, 3664299, 2557691, 1327176, 3011176, 3112561, 1004006, 1624584, 403534, 3602995, 1902346, 2901430, 1927354, 924025, 
4045597, 1654486, 1916225, 871954, 1180767, 3234430, 2409717, 3890788, 4197121, 151478, 801801, 3708044, 2567986, 3561609, 426783, 569208, 
1809453, 3236752, 2768026, 1151268, 443551, 301051, 3108584, 1554794, 1330423, 52367, 3991941, 3301299, 2349656, 1597758, 1438001, 632440, 
287830, 983850, 4174588, 676150, 48540, 2510140, 493846, 3999093, 3219668, 435598, 4006133, 3600162, 1099210, 1022256, 3829192, 3702123, 
2409970, 3985334, 1950414, 4104558, 2315120, 1361328, 3661545, 1092060, 3906808, 883401, 3625592, 1174067, 2287528, 4168793, 3265043, 1846113, 
1069857, 2839099, 2227672, 3548549, 2413133, 944126, 113633, 3408137, 248770, 2265785, 3137623, 1649085, 2804160, 3531167, 102525, 2665586, 
1618725, 2343229, 2395802, 2897351, 3314773, 3088390, 3178953, 1448224, 2257844, 42544, 2528686, 114965, 3412909, 2833541, 453332, 2609432, 
3033468, 2098783, 2817166, 759780, 3922382, 2014627, 320319, 1220152, 2038719, 3185328, 2021290, 4154200, 3945204, 400821, 18323, 3289152, 
92265, 3486946, 4092349, 4061764, 3125308, 2733357, 2631868, 2211356, 1909991, 3929410, 4096298, 416434, 2442426, 1100083, 3981631, 1560066, 
2479932, 434552, 3145437, 1326938, 907524, 2181622, 4029246, 831462, 2627803, 1257580, 2746597, 3235757, 3230939, 25302, 3534793, 4170750, 
};
