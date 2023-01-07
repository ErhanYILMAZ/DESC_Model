% This m file contains parameters for DESC model with 2-RC circuits.
% To estimate parameters UDDS data used in the data folder.

% Copyright (c) 2021 by Erhan YILMAZ (https://orcid.org/0000-0003-4788-9022)
% This work is licensed under CC BY-SA (Attribution-ShareAlike).
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.

%% DESC Model UDDS OCV Values 

OCV0 = [2.87330000000000,2.97390000000000,3.02730000000000,3.07211190476191,3.11221095238095,3.13925000000000,3.16411095238095,3.18520285714286,3.20393714285714,3.22150476190476,3.23921476190476,3.16550952380952,3.19831809523810,3.21685666666667,3.18136809523810,3.19385238095238,3.19364095238095,3.19869047619048,3.21232523809524,3.20610333333333,3.20839761904762,3.19478000000000,3.21222619047619,3.23085857142857,3.23440142857143,3.22434285714286,3.22882190476191,3.26420904761905,3.21001047619048,3.22845666666667,3.22814857142857,3.21842666666667,3.22039142857143,3.22764142857143,3.24078809523810,3.22095523809524,3.22720523809524,3.22444523809524,3.24785428571429,3.24101000000000,3.23869809523810,3.23828904761905,3.24480904761905,3.25043428571429,3.24555904761905,3.25469476190476,3.25571952380952,3.26782190476191,3.25244428571429,3.26076285714286,3.25984476190476,3.26954047619048,3.25929190476191,3.25982761904762,3.26049619047619,3.27900904761905,3.27122285714286,3.27357523809524,3.27076952380952,3.27574190476190,3.27916952380952,3.27445904761905,3.27374190476191,3.27549095238095,3.28017619047619,3.27491285714286,3.28400047619048,3.27851047619048,3.28660333333333,3.27957095238095,3.28119809523809,3.28336952380952,3.28908904761905,3.28391285714286,3.28370095238095,3.28187619047619,3.29113238095238,3.29122142857143,3.28957333333333,3.28899619047619,3.28996666666667,3.29493666666667,3.29284047619048,3.29312095238095,3.28849333333333,3.29318904761905,3.28700904761905,3.29632428571429,3.29365619047619,3.29959047619048,3.29727571428571,3.29714047619048,3.29796047619048,3.29988761904762,3.29840285714286,3.29440095238095,3.29180000000000,3.29540333333333,3.30613619047619,3.30330666666667,3.30419190476191,3.29914904761905,3.30455476190476,3.30398523809524,3.30311380952381,3.29542666666667,3.30154809523810,3.29670142857143,3.30510142857143,3.30628142857143,3.30824142857143,3.30833095238095,3.30425142857143,3.30633761904762,3.30606714285714,3.30991095238095,3.30620571428571,3.30346952380952,3.30036476190476,3.31436428571429,3.31146476190476,3.31593285714286,3.30040619047619,3.31185761904762,3.30414666666667,3.31443333333333,3.30319523809524,3.31365619047619,3.30586000000000,3.31326380952381,3.30735523809524,3.31405142857143,3.31326047619048,3.31594761904762,3.31582380952381,3.31772666666667,3.31946904761905,3.31607761904762,3.32066619047619,3.31832904761905,3.32425333333333,3.31754380952381,3.32810714285714,3.31987142857143,3.32816619047619,3.32015000000000,3.33039761904762,3.32298857142857,3.33055714285714,3.32494380952381,3.33174238095238,3.32889809523810,3.33273571428571,3.33123809523810,3.33106523809524,3.33425904761905,3.33397714285714,3.33774428571429,3.33295380952381,3.34073952380952,3.33336190476191,3.34104904761905,3.33332809523810,3.34372476190476,3.33267809523809,3.34138761904762,3.33890666666667,3.34578142857143,3.33936619047619,3.34240619047619,3.34263095238095,3.34326761904762,3.34378333333333,3.33974714285714,3.34960285714286,3.34529523809524,3.35212380952381,3.33998095238095,3.35214095238095,3.34614095238095,3.35656095238095,3.34539619047619,3.36236714285714,3.36014809523810,3.36144380952381,3.36187095238095,3.36240047619048,3.36304190476191,3.36348047619048,3.36399714285714,3.36448285714286,3.36495476190476,3.36544000000000,3.36573523809524,3.36616285714286,3.36592952380952,3.36676428571429,3.36111857142857,3.38028142857143,3.33446142857143,3.67780285714286];
OCVRel=[3.28069823170667e-19,1.34913183219078e-17,7.19857885040212e-18,-0.000481428571428580,-0.00133971428571428,-0.00161000000000002,-0.00193971428571430,-0.00166514285714289,-0.00158685714285715,-0.00166857142857142,-0.00182657142857143,-1.71428571428737e-05,-0.000212571428571443,-0.000181999999999999,0.000477428571428569,0.000165714285714264,0.000186285714285733,9.71428571428720e-05,-8.54285714285751e-05,-5.99999999996913e-06,-7.57142857143069e-05,0.000176000000000025,-0.000127142857142842,-0.000505428571428574,-0.000282571428571462,-0.000157142857142881,-0.000219428571428591,-0.000656285714285698,6.11428571428764e-05,-0.000142000000000015,-0.000147428571428544,5.20000000000138e-05,-4.45714285714262e-05,-3.45714285714177e-05,-5.85714285714356e-05,0.000520571428571415,0.000310571428571429,0.000438571428571433,-5.77142857143183e-05,0.000182000000000037,0.000243428571428593,0.000219714285714288,6.37142857142688e-05,-1.71428571428655e-06,0.000153714285714288,0.000189428571428581,0.000204857142857141,-5.94285714285823e-05,0.000260285714285737,0.000106857142857152,0.000179428571428571,-7.28571428571647e-05,0.000174571428571416,7.02857142857153e-05,0.000166857142857144,-7.62857142857085e-05,0.000238857142857155,0.000144571428571459,0.000214857142857137,0.000104571428571425,9.48571428571530e-05,0.000213714285714277,0.000184571428571405,0.000196285714285689,0.000142857142857154,0.000296857142857157,0.000259142857142869,0.000341142857142876,0.000173999999999982,0.000312285714285715,0.000283428571428595,0.000274857142857166,0.000119714285714279,0.000216857142857114,0.000218285714285690,0.000282857142857163,0.000221714285714280,0.000321428571428556,0.000267999999999985,0.000266857142857128,0.000260000000000016,0.000174000000000032,0.000207142857142843,0.000122285714285738,0.000212000000000004,0.000159714285714284,0.000223714285714306,0.000196285714285729,0.000318857142857152,0.000177142857142864,0.000203714285714295,0.000207142857142883,0.000191142857142873,0.000122285714285714,9.48571428571645e-05,0.000158285714285754,0.000179999999999969,0.000134000000000007,0.000154857142857124,0.000188000000000011,0.000134571428571414,0.000231714285714296,0.000121428571428564,0.000126571428571424,7.51428571428628e-05,0.000211999999999993,0.000113428571428573,0.000117428571428570,7.74285714285804e-05,0.000213428571428586,0.000145428571428577,0.000104285714285742,0.000167428571428586,0.000112285714285736,0.000139142857142844,2.02857142856830e-05,8.97142857143165e-05,7.48571428571617e-05,0.000123428571428566,8.42857142857417e-05,0.000163428571428597,8.57142857132929e-07,0.000388857142857112,0.000136285714285698,0.000436000000000012,8.00000000000102e-05,0.000388571428571394,0.000138857142857117,0.000332000000000002,0.000325142857142875,0.000480571428571411,0.000327428571428539,0.000371142857142840,0.000274285714285722,0.000317142857142870,0.000271999999999982,0.000215714285714314,0.000340285714285695,0.000200857142857140,0.000447714285714300,0.000424000000000044,0.000501142857142851,0.000207142857142876,0.000411428571428572,0.000220857142857156,0.000470000000000016,8.42857142857369e-05,0.000260571428571419,0.000177142857142869,0.000261142857142866,0.000263714285714285,0.000403428571428585,0.000275714285714269,0.000311428571428552,0.000282571428571437,0.000193714285714291,0.000201142857142878,4.02857142857207e-05,0.000163142857142829,-7.11428571428548e-05,6.85714285714191e-05,0.000171714285714285,0.000429428571428563,7.54285714285783e-05,0.000339428571428585,0.000122285714285751,0.000248000000000016,-2.65714285714159e-05,6.08571428571378e-05,6.88571428571347e-05,-9.57142857142537e-05,-4.17142857142528e-05,0.000209999999999983,0.000295142857142856,5.48571428571780e-05,0.000128571428571416,-8.28571428571514e-05,0.000294285714285723,-0.000113714285714299,6.28571428571398e-06,-0.000369714285714302,-0.000153142857142867,-0.000240857142857145,-0.000166571428571451,-0.000198857142857130,-0.000207714285714279,-0.000220857142857153,-0.000235428571428595,-0.000244857142857142,-0.000254857142857141,-0.000269142857142889,-0.000278571428571404,-0.000292000000000005,-0.000303428571428560,-0.000313142857142855,-0.000313142857142881,-0.000315714285714285,-0.000173428571428577,-0.000226571428571405,0.00120942857142861,-0.00116514285714285];

%% DESC Model UDDS -5° C

T=-5;
M = 0.053406;
R0cha = 0.0088025;
R0dch = 0.0057418;
R1 = 0.0014174;
R2 = 0.0071788;
RC1 = 6.1635;
RC2 = 255.97;
g = 80.216;
Q=  109.0828;

%% DESC Model UDDS 5° C

T=5;
M = 0.042486;
R0cha = 0.0043191;
R0dch = 0.0037063;
R1 = 0.0012043;
R2 = 0.0051726;
RC1 = 9.2128;
RC2 = 284.7;
g = 121.02;
Q = 109.1659;
    
%% DESC Model UDDS 15° C

T=15;
M = 0.029216;
R0cha = 0.0021067;
R0dch = 0.0019291;
R1 = 0.00091722;
R2 = 0.0026572;
RC1 = 17.347;
RC2 = 336.76;
g = 33.94;
Q=109.3002;

%% DESC Model UDDS 25° C
        
T=25;
M = 0.024564;
R0cha = 0.0012828;
R0dch = 0.0012426;
R1 = 0.00063893;
R2 = 0.0019287;
RC1 = 16.117;
RC2 = 331.3;
g = 31.283;
Q=109.3742;

%% DESC Model UDDS 35° C
         
T=35;
M = 0.021225;
R0cha = 0.00088795;
R0dch = 0.0009058;
R1 = 0.00039402;
R2 = 0.001625;
RC1 = 12.696;
RC2 = 307.16;
g = 126.45;
Q=109.5282;

%% DESC Model UDDS 45° C

T=45;
M = 0.022885;
R0cha = 0.00071734;
R0dch = 0.00071758;
R1 = 0.00039333;
R2 = 0.001206;
RC1 = 14.822;
RC2 = 399.55;
g = 54.84;
Q=109.6644;

%% Identified OCV values for test temperature values
% These values used to get OCV0 and OCVRel valuese with least square regression

OCV_N05 = [2.87330000000000,2.97390000000000,3.02730000000000,3.06810000000000,3.10170000000000,3.12990000000000,3.15460000000000,3.17660000000000,3.19590000000000,3.21300000000000,3.23400000000000,3.14950000000000,3.18800000000000,3.21130000000000,3.17100000000000,3.18420000000000,3.18410000000000,3.18930000000000,3.20530000000000,3.19790000000000,3.20000000000000,3.18410000000000,3.20360000000000,3.22420000000000,3.22960000000000,3.21730000000000,3.22200000000000,3.21830000000000,3.22660000000000,3.23200000000000,3.22970000000000,3.22590000000000,3.22590000000000,3.23390000000000,3.25290000000000,3.24820000000000,3.25080000000000,3.25030000000000,3.24810000000000,3.25910000000000,3.25780000000000,3.25660000000000,3.24940000000000,3.25460000000000,3.26490000000000,3.27830000000000,3.26930000000000,3.27310000000000,3.26960000000000,3.27170000000000,3.27840000000000,3.27460000000000,3.27150000000000,3.26710000000000,3.27160000000000,3.28380000000000,3.28420000000000,3.28460000000000,3.28330000000000,3.27890000000000,3.28570000000000,3.28620000000000,3.28230000000000,3.27750000000000,3.27810000000000,3.28190000000000,3.29280000000000,3.28750000000000,3.28940000000000,3.28640000000000,3.28590000000000,3.29340000000000,3.29000000000000,3.28750000000000,3.28350000000000,3.28470000000000,3.29420000000000,3.29660000000000,3.29640000000000,3.29550000000000,3.29170000000000,3.29570000000000,3.29870000000000,3.29430000000000,3.29030000000000,3.28990000000000,3.28580000000000,3.30180000000000,3.30260000000000,3.30240000000000,3.30040000000000,3.29740000000000,3.30430000000000,3.30150000000000,3.29820000000000,3.29370000000000,3.29060000000000,3.29560000000000,3.30900000000000,3.30770000000000,3.30730000000000,3.30300000000000,3.30420000000000,3.30950000000000,3.30450000000000,3.30010000000000,3.29890000000000,3.28860000000000,3.30370000000000,3.31340000000000,3.31150000000000,3.30980000000000,3.30470000000000,3.31200000000000,3.31120000000000,3.30690000000000,3.30340000000000,3.29720000000000,3.29710000000000,3.31540000000000,3.31520000000000,3.32280000000000,3.30070000000000,3.30980000000000,3.29190000000000,3.31200000000000,3.29610000000000,3.31260000000000,3.29570000000000,3.31290000000000,3.30100000000000,3.31260000000000,3.30340000000000,3.31140000000000,3.31010000000000,3.31570000000000,3.31240000000000,3.30950000000000,3.31940000000000,3.31550000000000,3.32140000000000,3.30780000000000,3.32630000000000,3.31410000000000,3.32590000000000,3.30970000000000,3.32880000000000,3.31600000000000,3.32700000000000,3.31610000000000,3.32820000000000,3.32340000000000,3.32960000000000,3.32390000000000,3.32660000000000,3.33250000000000,3.33120000000000,3.33250000000000,3.32500000000000,3.33940000000000,3.33020000000000,3.33840000000000,3.32510000000000,3.34360000000000,3.32940000000000,3.34000000000000,3.32930000000000,3.34350000000000,3.33490000000000,3.33970000000000,3.33660000000000,3.34200000000000,3.34340000000000,3.34000000000000,3.34480000000000,3.34040000000000,3.35220000000000,3.33740000000000,3.35060000000000,3.34200000000000,3.36070000000000,3.35160000000000,3.36310000000000,3.34760000000000,3.35260000000000,3.35230000000000,3.35320000000000,3.35400000000000,3.35460000000000,3.35530000000000,3.35590000000000,3.35670000000000,3.35740000000000,3.35790000000000,3.35860000000000,3.35870000000000,3.36040000000000,3.35460000000000,3.37420000000000,3.31050000000000,3.67390000000000];
OCV_P05 = [2.87330000000000,2.97390000000000,3.02730000000000,3.06810000000000,3.10170000000000,3.12990000000000,3.15460000000000,3.17660000000000,3.19590000000000,3.21300000000000,3.23400000000000,3.14950000000000,3.18800000000000,3.21130000000000,3.17100000000000,3.18420000000000,3.18410000000000,3.18930000000000,3.20530000000000,3.19790000000000,3.20000000000000,3.18410000000000,3.20360000000000,3.22420000000000,3.22960000000000,3.21730000000000,3.22200000000000,3.21830000000000,3.22660000000000,3.23200000000000,3.22970000000000,3.22590000000000,3.22590000000000,3.23610000000000,3.23750000000000,3.19180000000000,3.21050000000000,3.19650000000000,3.24950000000000,3.22430000000000,3.21800000000000,3.22140000000000,3.24430000000000,3.25210000000000,3.22380000000000,3.24030000000000,3.23770000000000,3.26350000000000,3.23270000000000,3.25240000000000,3.24050000000000,3.26770000000000,3.24530000000000,3.25930000000000,3.24900000000000,3.27580000000000,3.25960000000000,3.26300000000000,3.26150000000000,3.27620000000000,3.27530000000000,3.26430000000000,3.27230000000000,3.27270000000000,3.28380000000000,3.26730000000000,3.27980000000000,3.27110000000000,3.28630000000000,3.27320000000000,3.28330000000000,3.27590000000000,3.28900000000000,3.28040000000000,3.28530000000000,3.28260000000000,3.29110000000000,3.28880000000000,3.28600000000000,3.28980000000000,3.28950000000000,3.29540000000000,3.28650000000000,3.29500000000000,3.28740000000000,3.29750000000000,3.28940000000000,3.29740000000000,3.28950000000000,3.29730000000000,3.29430000000000,3.29830000000000,3.29430000000000,3.29810000000000,3.29950000000000,3.29850000000000,3.29980000000000,3.29660000000000,3.30390000000000,3.29840000000000,3.30420000000000,3.29520000000000,3.30470000000000,3.30000000000000,3.30680000000000,3.29600000000000,3.30370000000000,3.30390000000000,3.30740000000000,3.30100000000000,3.30320000000000,3.30820000000000,3.30730000000000,3.30720000000000,3.30220000000000,3.31210000000000,3.30750000000000,3.31200000000000,3.30170000000000,3.31300000000000,3.30890000000000,3.31560000000000,3.30360000000000,3.31250000000000,3.31250000000000,3.31730000000000,3.30980000000000,3.31310000000000,3.31750000000000,3.31870000000000,3.31720000000000,3.31580000000000,3.32190000000000,3.32030000000000,3.32430000000000,3.31850000000000,3.32590000000000,3.32270000000000,3.32920000000000,3.32300000000000,3.32820000000000,3.32590000000000,3.33290000000000,3.32780000000000,3.33090000000000,3.32980000000000,3.33490000000000,3.33320000000000,3.33380000000000,3.33420000000000,3.33660000000000,3.33790000000000,3.33630000000000,3.33870000000000,3.33780000000000,3.34260000000000,3.33720000000000,3.34240000000000,3.33950000000000,3.34540000000000,3.33830000000000,3.34340000000000,3.34220000000000,3.34770000000000,3.34080000000000,3.34330000000000,3.34650000000000,3.34840000000000,3.34560000000000,3.34350000000000,3.35100000000000,3.34880000000000,3.35150000000000,3.34150000000000,3.35450000000000,3.34970000000000,3.35580000000000,3.34000000000000,3.35440000000000,3.35400000000000,3.36230000000000,3.35030000000000,3.36780000000000,3.38340000000000,3.37990000000000,3.38150000000000,3.38170000000000,3.38220000000000,3.38270000000000,3.38310000000000,3.38360000000000,3.38380000000000,3.38420000000000,3.38430000000000,3.38450000000000,3.38390000000000,3.38380000000000,3.37810000000000,3.39260000000000,3.35970000000000,3.65820000000000];
OCV_P15 = [2.87330000000000,2.97390000000000,3.02730000000000,3.06810000000000,3.10170000000000,3.12990000000000,3.15460000000000,3.17660000000000,3.19590000000000,3.21300000000000,3.22710000000000,3.23600000000000,3.23860000000000,3.23870000000000,3.23960000000000,3.24000000000000,3.24060000000000,3.24100000000000,3.24180000000000,3.24290000000000,3.24390000000000,3.24550000000000,3.24700000000000,3.24870000000000,3.25050000000000,3.25310000000000,3.25540000000000,3.42820000000000,3.15440000000000,3.21840000000000,3.22540000000000,3.19770000000000,3.20480000000000,3.20510000000000,3.22070000000000,3.21570000000000,3.20730000000000,3.21870000000000,3.24790000000000,3.23070000000000,3.23340000000000,3.22880000000000,3.24080000000000,3.24230000000000,3.23830000000000,3.22740000000000,3.25750000000000,3.26380000000000,3.24970000000000,3.25180000000000,3.25100000000000,3.26180000000000,3.25730000000000,3.24840000000000,3.25660000000000,3.27310000000000,3.26510000000000,3.26770000000000,3.26160000000000,3.27120000000000,3.27240000000000,3.26720000000000,3.26040000000000,3.27920000000000,3.28260000000000,3.27580000000000,3.27510000000000,3.27460000000000,3.28360000000000,3.27910000000000,3.27220000000000,3.27600000000000,3.28900000000000,3.28410000000000,3.28580000000000,3.27920000000000,3.28860000000000,3.28840000000000,3.28310000000000,3.27750000000000,3.29090000000000,3.29470000000000,3.29150000000000,3.28890000000000,3.28880000000000,3.29670000000000,3.29030000000000,3.28530000000000,3.28360000000000,3.29840000000000,3.29730000000000,3.29720000000000,3.29100000000000,3.29980000000000,3.29780000000000,3.29180000000000,3.28570000000000,3.29530000000000,3.30410000000000,3.30210000000000,3.29820000000000,3.29890000000000,3.30550000000000,3.29790000000000,3.29470000000000,3.28640000000000,3.30530000000000,3.30670000000000,3.30570000000000,3.29900000000000,3.30870000000000,3.30530000000000,3.29970000000000,3.29220000000000,3.30040000000000,3.31210000000000,3.31030000000000,3.30550000000000,3.30760000000000,3.31300000000000,3.30680000000000,3.30140000000000,3.29560000000000,3.31500000000000,3.31500000000000,3.31400000000000,3.30850000000000,3.31850000000000,3.31480000000000,3.31170000000000,3.30540000000000,3.31480000000000,3.32120000000000,3.32110000000000,3.31560000000000,3.32030000000000,3.32450000000000,3.32080000000000,3.31380000000000,3.32130000000000,3.32830000000000,3.32640000000000,3.32570000000000,3.32220000000000,3.33140000000000,3.32740000000000,3.32650000000000,3.32210000000000,3.33310000000000,3.33210000000000,3.33390000000000,3.32770000000000,3.33400000000000,3.33530000000000,3.33150000000000,3.32560000000000,3.33370000000000,3.34040000000000,3.33890000000000,3.33730000000000,3.33460000000000,3.34290000000000,3.33600000000000,3.33670000000000,3.32680000000000,3.34130000000000,3.34500000000000,3.34440000000000,3.33820000000000,3.34510000000000,3.34390000000000,3.33900000000000,3.33140000000000,3.33560000000000,3.35020000000000,3.34830000000000,3.34520000000000,3.34350000000000,3.35060000000000,3.34290000000000,3.34050000000000,3.32890000000000,3.35250000000000,3.34960000000000,3.35080000000000,3.35080000000000,3.35100000000000,3.35150000000000,3.35160000000000,3.35180000000000,3.35200000000000,3.35220000000000,3.35240000000000,3.35250000000000,3.35270000000000,3.35250000000000,3.35320000000000,3.34940000000000,3.37240000000000,3.35170000000000,3.68930000000000];
OCV_P25 = [2.87330000000000,2.97390000000000,3.02730000000000,3.06810000000000,3.10170000000000,3.12990000000000,3.15460000000000,3.17660000000000,3.19590000000000,3.21300000000000,3.21230000000000,3.13620000000000,3.17800000000000,3.20270000000000,3.16390000000000,3.17690000000000,3.17620000000000,3.18350000000000,3.19760000000000,3.18970000000000,3.19200000000000,3.17610000000000,3.19610000000000,3.21460000000000,3.22000000000000,3.20780000000000,3.21160000000000,3.20790000000000,3.21720000000000,3.22090000000000,3.21910000000000,3.21590000000000,3.21640000000000,3.22570000000000,3.24280000000000,3.23940000000000,3.24180000000000,3.24100000000000,3.24090000000000,3.25080000000000,3.24940000000000,3.24850000000000,3.24190000000000,3.24810000000000,3.25730000000000,3.27030000000000,3.26270000000000,3.26580000000000,3.26270000000000,3.26660000000000,3.27150000000000,3.26850000000000,3.26600000000000,3.26180000000000,3.26790000000000,3.27830000000000,3.28020000000000,3.28070000000000,3.27940000000000,3.27720000000000,3.28400000000000,3.28450000000000,3.28220000000000,3.27870000000000,3.28100000000000,3.28420000000000,3.29420000000000,3.29050000000000,3.29190000000000,3.28890000000000,3.29010000000000,3.29560000000000,3.29260000000000,3.29100000000000,3.28700000000000,3.28940000000000,3.29690000000000,3.30050000000000,3.30010000000000,3.29890000000000,3.29600000000000,3.30010000000000,3.30190000000000,3.29820000000000,3.29450000000000,3.29490000000000,3.29060000000000,3.30490000000000,3.30710000000000,3.30630000000000,3.30380000000000,3.30270000000000,3.30770000000000,3.30510000000000,3.30210000000000,3.29810000000000,3.29580000000000,3.29900000000000,3.31290000000000,3.31160000000000,3.31070000000000,3.30720000000000,3.30910000000000,3.31250000000000,3.30810000000000,3.30430000000000,3.30370000000000,3.29400000000000,3.30700000000000,3.31810000000000,3.31560000000000,3.31360000000000,3.31050000000000,3.31570000000000,3.31520000000000,3.31130000000000,3.30850000000000,3.30310000000000,3.30150000000000,3.32020000000000,3.32050000000000,3.32010000000000,3.31660000000000,3.31710000000000,3.32200000000000,3.31940000000000,3.31690000000000,3.31670000000000,3.31010000000000,3.32010000000000,3.32720000000000,3.32720000000000,3.32560000000000,3.32300000000000,3.32780000000000,3.32960000000000,3.32690000000000,3.32710000000000,3.32540000000000,3.33100000000000,3.33600000000000,3.33430000000000,3.33540000000000,3.33270000000000,3.33310000000000,3.33750000000000,3.33510000000000,3.33270000000000,3.33690000000000,3.33000000000000,3.33940000000000,3.34350000000000,3.34330000000000,3.34220000000000,3.33990000000000,3.34280000000000,3.34380000000000,3.33990000000000,3.33820000000000,3.33680000000000,3.33300000000000,3.34780000000000,3.35060000000000,3.34920000000000,3.34740000000000,3.34600000000000,3.35020000000000,3.34740000000000,3.34250000000000,3.34590000000000,3.33560000000000,3.34020000000000,3.35580000000000,3.35390000000000,3.35330000000000,3.34970000000000,3.35100000000000,3.35420000000000,3.34910000000000,3.34560000000000,3.34190000000000,3.33740000000000,3.35400000000000,3.35370000000000,3.35410000000000,3.35420000000000,3.35440000000000,3.35450000000000,3.35450000000000,3.35480000000000,3.35480000000000,3.35500000000000,3.35500000000000,3.35500000000000,3.35520000000000,3.35490000000000,3.35550000000000,3.35390000000000,3.37280000000000,3.37180000000000,3.64910000000000];
OCV_P35 = [2.87330000000000,2.97390000000000,3.02730000000000,3.06810000000000,3.09190000000000,3.06890000000000,3.05480000000000,3.09750000000000,3.11810000000000,3.13200000000000,3.14800000000000,3.15460000000000,3.18440000000000,3.20690000000000,3.19460000000000,3.19690000000000,3.19680000000000,3.19880000000000,3.20590000000000,3.20160000000000,3.20080000000000,3.19810000000000,3.20190000000000,3.20660000000000,3.22280000000000,3.21490000000000,3.21750000000000,3.21600000000000,3.22080000000000,3.22410000000000,3.22210000000000,3.22490000000000,3.22000000000000,3.22970000000000,3.24230000000000,3.24530000000000,3.24540000000000,3.24590000000000,3.24660000000000,3.25150000000000,3.25090000000000,3.24980000000000,3.24970000000000,3.25180000000000,3.25340000000000,3.26860000000000,3.26620000000000,3.26740000000000,3.26590000000000,3.26810000000000,3.27010000000000,3.26780000000000,3.26860000000000,3.26470000000000,3.26870000000000,3.27750000000000,3.28310000000000,3.28220000000000,3.28230000000000,3.28190000000000,3.28450000000000,3.28480000000000,3.28290000000000,3.28360000000000,3.28500000000000,3.28630000000000,3.29550000000000,3.29360000000000,3.29400000000000,3.29270000000000,3.29390000000000,3.29510000000000,3.29310000000000,3.29240000000000,3.29270000000000,3.29320000000000,3.29970000000000,3.30410000000000,3.30160000000000,3.30130000000000,3.30080000000000,3.30140000000000,3.30070000000000,3.29750000000000,3.29740000000000,3.29840000000000,3.29440000000000,3.30510000000000,3.30870000000000,3.30700000000000,3.30590000000000,3.30580000000000,3.30620000000000,3.30410000000000,3.30130000000000,3.30170000000000,3.29860000000000,3.29950000000000,3.31270000000000,3.31150000000000,3.31040000000000,3.30930000000000,3.30900000000000,3.30940000000000,3.30610000000000,3.30560000000000,3.30480000000000,3.29850000000000,3.30700000000000,3.31670000000000,3.31420000000000,3.31310000000000,3.31190000000000,3.31230000000000,3.31200000000000,3.30940000000000,3.30940000000000,3.30470000000000,3.30280000000000,3.31760000000000,3.31940000000000,3.31830000000000,3.31680000000000,3.31670000000000,3.31820000000000,3.31650000000000,3.31750000000000,3.31790000000000,3.31630000000000,3.32500000000000,3.32600000000000,3.32670000000000,3.32590000000000,3.32600000000000,3.32690000000000,3.32730000000000,3.32610000000000,3.32880000000000,3.32910000000000,3.33550000000000,3.34010000000000,3.33540000000000,3.33650000000000,3.33580000000000,3.33660000000000,3.33610000000000,3.33320000000000,3.33290000000000,3.33720000000000,3.33320000000000,3.34120000000000,3.34520000000000,3.34370000000000,3.34290000000000,3.34280000000000,3.34280000000000,3.34130000000000,3.33810000000000,3.33860000000000,3.33810000000000,3.33500000000000,3.34730000000000,3.35010000000000,3.34800000000000,3.34720000000000,3.34700000000000,3.34670000000000,3.34430000000000,3.34160000000000,3.34510000000000,3.33760000000000,3.34120000000000,3.35420000000000,3.35250000000000,3.35170000000000,3.35040000000000,3.34970000000000,3.35060000000000,3.34730000000000,3.34700000000000,3.34420000000000,3.34110000000000,3.35320000000000,3.35360000000000,3.35360000000000,3.35380000000000,3.35380000000000,3.35390000000000,3.35400000000000,3.35420000000000,3.35410000000000,3.35420000000000,3.35410000000000,3.35390000000000,3.35380000000000,3.35340000000000,3.35420000000000,3.35370000000000,3.37170000000000,3.37810000000000,3.65300000000000];
OCV_P45 = [2.87330000000000,2.97390000000000,3.02730000000000,3.03440000000000,3.01380000000000,3.05380000000000,3.07870000000000,3.10750000000000,3.13150000000000,3.14480000000000,3.16070000000000,3.16520000000000,3.18740000000000,3.20840000000000,3.20540000000000,3.20080000000000,3.20240000000000,3.20190000000000,3.20780000000000,3.20590000000000,3.20460000000000,3.20190000000000,3.20590000000000,3.20620000000000,3.22000000000000,3.21680000000000,3.21810000000000,3.21780000000000,3.22180000000000,3.22630000000000,3.22520000000000,3.22650000000000,3.22400000000000,3.23120000000000,3.24150000000000,3.24780000000000,3.24470000000000,3.24690000000000,3.24720000000000,3.25150000000000,3.25190000000000,3.25100000000000,3.25040000000000,3.25350000000000,3.25410000000000,3.26600000000000,3.26550000000000,3.26620000000000,3.26530000000000,3.26680000000000,3.26910000000000,3.26810000000000,3.26800000000000,3.26610000000000,3.26920000000000,3.27640000000000,3.28380000000000,3.28060000000000,3.28230000000000,3.28160000000000,3.28450000000000,3.28540000000000,3.28450000000000,3.28480000000000,3.28770000000000,3.28960000000000,3.29770000000000,3.29470000000000,3.29530000000000,3.29460000000000,3.29580000000000,3.29720000000000,3.29520000000000,3.29410000000000,3.29410000000000,3.29610000000000,3.30290000000000,3.30750000000000,3.30240000000000,3.30300000000000,3.30210000000000,3.30320000000000,3.30260000000000,3.29950000000000,3.29800000000000,3.30090000000000,3.29840000000000,3.30700000000000,3.30870000000000,3.30740000000000,3.30640000000000,3.30630000000000,3.30720000000000,3.30540000000000,3.30290000000000,3.30160000000000,3.30190000000000,3.30250000000000,3.31280000000000,3.31110000000000,3.31050000000000,3.30910000000000,3.30940000000000,3.30980000000000,3.30750000000000,3.30560000000000,3.30650000000000,3.30260000000000,3.30910000000000,3.31510000000000,3.31370000000000,3.31250000000000,3.31150000000000,3.31210000000000,3.31210000000000,3.31010000000000,3.30890000000000,3.30730000000000,3.30630000000000,3.31710000000000,3.31760000000000,3.31750000000000,3.31580000000000,3.31640000000000,3.31760000000000,3.31700000000000,3.31700000000000,3.31980000000000,3.32060000000000,3.33020000000000,3.32500000000000,3.32650000000000,3.32610000000000,3.32680000000000,3.32830000000000,3.32760000000000,3.32690000000000,3.32840000000000,3.33120000000000,3.33740000000000,3.34240000000000,3.33560000000000,3.33670000000000,3.33600000000000,3.33760000000000,3.33680000000000,3.33400000000000,3.33230000000000,3.33660000000000,3.33540000000000,3.34280000000000,3.34410000000000,3.34260000000000,3.34180000000000,3.34170000000000,3.34250000000000,3.34080000000000,3.33800000000000,3.33710000000000,3.33890000000000,3.33730000000000,3.34710000000000,3.34750000000000,3.34620000000000,3.34520000000000,3.34540000000000,3.34550000000000,3.34350000000000,3.34070000000000,3.34340000000000,3.33960000000000,3.34340000000000,3.35160000000000,3.35040000000000,3.34970000000000,3.34870000000000,3.34890000000000,3.34950000000000,3.34720000000000,3.34610000000000,3.34540000000000,3.34470000000000,3.35470000000000,3.35300000000000,3.35380000000000,3.35370000000000,3.35380000000000,3.35390000000000,3.35410000000000,3.35420000000000,3.35420000000000,3.35440000000000,3.35450000000000,3.35440000000000,3.35460000000000,3.35460000000000,3.35560000000000,3.35620000000000,3.37080000000000,3.38010000000000,3.60350000000000];

    