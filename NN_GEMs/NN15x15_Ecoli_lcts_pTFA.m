function [flux,flux_bio] = NN15x15_Ecoli_lcts_o2_meth_ac_gal_20K_pTFA(vshrinkage,x,Ncells,num_metInert,~,~,~)

%x=[lcts,o2,meth,ac,gal]
%input for NN: x_cut=[lcts,o2,meth,ac,gal]
%Value below which an uptake flux [mmol gDW^-1 h^-1] is negligible because
%the biomass predicted by TFA is zero.
% negl_l=10^-3;    %Lactose
% negl_o=0;   %Oxygen
% negl_m=10^-5;    %Methionine
% negl_a=10^-3;    %Acetate
% negl_g=10^-3;    %Galactose
x_cut = [x(:,1).*(x(:,1)<-10^-3),x(:,2),x(:,3).*(x(:,3)<-10^-5),x(:,4).*(x(:,4)<-10^-3),x(:,5).*(x(:,5)<-10^-3)];

%E. coli is not divided in subpopulations 
x_sqrt = sqrt([-x_cut]');


% ===== NEURAL NETWORK CONSTANTS =====

% Layer 1
b1 = [-0.32607280319851439376;0.90132493644465050409;0.70586781739483595022;2.0506921654889205264;2.1652198189565781128;-0.33474989451958203635;-1.3633584406426022539;-0.3696099460104174117;0.56916754787397627613;0.77852397984137056142;0.67241797453432317067;-0.2450738769637219594;-0.22033774180062112857;-1.6429748099140002449;0.92144274112748025818];
IW1_1 = [0.37849178799315913446 0.062475677377128183143 -4.2767083244861705538 -0.0012672673447689151126 -0.057931275142581374615;-0.041301758178887919282 -0.93806147971193876156 8.4790335053445478053 0.0035252551907242663833 -0.027784122495657453955;-0.053588789367684415421 -0.95348545088680214388 8.7464675170654473391 0.0034604872875669424187 -0.035276771055136213473;1.4854748493481626515 -1.1291834954244952538 -2.2253200719411676545 0.095821055166951399351 0.92442713013016730716;-2.6412642496873868225 -2.6695677446949206413 41.634665963722923721 0.010600460644199128854 -0.17560430032197099637;0.36291735416110332269 0.056226160009158435149 -4.5582394246519788439 0.0010962663117655133994 0.016446193291457012614;0.946702071212811469 0.084452572627271751959 -5.1809127079190906429 -0.0086968627155679326807 0.085634839236688062303;-0.21200891041799141412 -0.020499461341386186752 1.2981232980005503741 -0.00054273531758322894787 -0.090895581725277785345;-0.027491480265920187342 0.18193212515196147128 -0.6846325134498828513 0.0015680070992486968857 -0.0011714328849612583826;0.080866829727019029495 -0.64632691212709081352 2.5300920933996988715 0.002055159112928065749 -0.01672059705284266104;0.16952862858626663312 0.046283053375967393717 -1.9467604383198939821 -0.00010617677541185201802 -0.0043735599191148987266;0.0054784150542180080554 -0.0011179669769464390793 0.04750378330764585405 -0.026889456800475813186 0.0022579040930282875781;0.0083612459940824446214 0.0020551208344947776568 0.7683050730450206034 0.02061676389514309829 -0.0054906195236764316139;-0.25547670042860937123 -0.0019260339511782717352 -0.20640150932122236505 -0.0018751155290203614864 0.13634272681651660175;-0.31597880482126705637 -0.041901661261364825384 3.0647188361237356347 0.0093556181620264352444 -0.4268382544670415446];

% Layer 2
b2 = [22.926111807455392011;-10.700228533546471965;-12.785060052099368022;1.1005044096409588583;145.72128821460015047;0.45646026254626592689;-11.580744493789165261;0.48859613582803279597;1.2155715136121281184;13.664588665388640365;-0.32801936307391277303;-168.75518029568564771;4.5572726893597970133;-18.76826144575347044;-60.459655282418182765];
LW2_1 = [7.9824280607437243518 0.078819359634947336324 0.21181886842972749707 -0.11929423121510451811 0.033757451464400303753 -6.9738018484363681893 -0.75809567299870939117 -1.3029094820878350802 -1.3640165760778890647 -0.10767933761158765549 5.3277694114891502863 1.9280246933624491934 7.5372319599370420917 26.155286710845285114 -0.9396516083475034975;0.78834913378831295017 16.095185084365503769 -11.915497787687192144 -0.054631158923748361633 0.20710820540861613059 -2.4731174377638627249 0.34825708750388689738 -0.71677119723836724141 0.30269786713211183082 0.13365192624004501587 -0.54693062590676799761 -1.227655679368826247 -1.2873265552718784477 -4.431706586936495107 -0.080774063706690130138;26.508656369054360624 -26.705576975806035733 22.775628756524206153 -1.4425573170625989317 0.4822914307973794612 -37.241491538192342148 2.6566898905253291119 -16.267105809477072853 -3.6628270033721870291 -5.9331912254106731197 6.6006668782819390984 -39.322804428884801098 -42.829012866059301246 4.2111893695149493766 -3.0561854182795404888;0.24775694219270788521 -0.25703230005142707615 0.25118742804324339923 -0.28306355557627549402 0.13121368901188557032 -0.61871060203127292709 0.24592197690367864626 -0.36713209350816145005 -0.84981396227552907785 -0.073336131985866512273 0.20929120157051300799 0.075022240893694400188 -0.021630787973858783824 0.85636326212149516124 0.021178792419258109236;65.465552574943856712 9.5454292569422278802 -7.1872387891123841186 -0.79884274819937794998 0.051383064007625996183 -80.458389065810649754 30.768325149317735878 56.225998614383144059 -22.27005486414945068 -3.1429869816801789817 65.773380103698684707 -28.46097842144234491 25.763457312690952961 115.83647288836628775 -37.539672111658838105;-0.031339090938655625451 0.094736693637075042318 -0.076778022914219881678 0.062224081204371957388 -0.033355949759835530621 0.072180507024057816157 -0.1030695858119410474 -0.0031301873524941003016 0.22055787481524011251 0.0091063670432326611348 -0.030117643572723694356 4.0636805384884384651 -0.10779384850215363578 -0.41123448688002306461 0.01031733753832824764;0.45699976250486407103 16.934053959074905293 -12.569948023768366951 -0.037460112780389298592 0.23173166581614354498 -2.2102638320459000454 -0.092565459569237784621 -0.52248775852052098401 0.24282949036653411756 0.085343109415617235403 -0.046830679417625399596 -0.80883649592612771873 -0.92785143794184288257 -5.2639453628570453958 -0.016996629926163586005;0.066523748978646352059 0.15930126929060448027 -0.11850605442899417308 0.1549937631137013605 -0.094930882452214399025 0.2185135609185561989 0.43920509935101231802 1.0528969157683580882 0.48454116425754967246 0.097336896893239016615 0.24002029143777947606 -0.19359795188016515333 -0.0095694117148158643754 -0.043711189360326647979 -0.081816372327213796045;0.23136549098492620313 0.090965623619763888441 -0.072843210479062558971 0.097141300211623787808 -0.049928112373234488641 -0.10680303379709721323 -0.15887061916059008748 -0.0028579715706216202165 0.21492640941160764956 0.022426582957565294879 0.11351827144563436667 4.0473710511070324358 0.39469167266924459758 0.4713883660549793575 -0.029070094523662468738;-26.983428873419661898 25.583548942038369489 -21.83637023359908369 1.5353241120182374857 -0.48854494513620877472 37.894568624196310225 -2.6815408078396449909 16.637925384035195719 3.370379154341839989 5.9094718280398250698 -6.4914642152758901261 39.800686770069319209 43.384308595947558729 -4.0079061063471339565 3.0879279480409937975;-0.037279490960289429768 -0.21813406904151763643 0.17550522103068169621 -0.20854202707542687789 0.1201079367571012968 -0.33825830071038381908 -0.38106089910781187902 -1.0863777100991316438 -0.6127414021151972312 -0.11176842662096199421 -0.28938504539647114111 -0.038217700598388082422 -0.24741298233802314077 0.19830169260827881783 0.088111526317498578553;-14.059540539471136356 -106.83589014625233915 89.539058692654734273 -4.1101014263519735081 2.7356957091504070512 -54.384733780815373905 -2.6139401798001524746 -31.325708267281353159 -66.271432982268606793 -38.851885968169185048 321.99034935361169119 -51.282175587096553215 -56.372149647963048835 4.7631299176364585435 1.3163847915380952891;1.5589425095781990205 -0.25874728458515561247 0.19697059042286574337 0.23525715135883190632 -0.10783992210051905614 -1.0086952311522869863 -0.35511242803468229257 0.050530987932679426178 -0.41369818625984577931 0.12148745328215163064 0.82132970321886578091 -1.513106832278375169 2.644769523833714775 5.3383256473727573876 -0.23968450353391057561;1.8052272593828815328 31.977286881666152851 -23.601971724279938059 -0.092506418971662354722 -0.48701329949034044553 -4.8300995429796813596 0.74902516656346418422 -1.4913658729622658594 0.70310670972322619754 0.2967118516948246687 -1.4111803952923145378 -2.5264254861655102147 -2.6071761437551130847 -7.8293358796428504931 -0.11522732480101958996;-16.634424945298622589 17.556514478962871095 -12.952278328461154899 0.051431135709839974346 -0.72929879930954644163 15.285224701267905445 1.2733701570786979396 4.6061719155110241175 4.1725850176826222437 0.35734723977176746645 -14.522234509626317234 -10.201125002058180513 -13.429752722546522747 -60.711473325497976816 2.1218220020967959449];

% Layer 3
b3 = [0.3111390394728638098;1.7059303151423084177;0.25769290808224132316;0.68796358714795224909;-12.508594674744443864;17.346603004199131703;-13.029121904207194405;0.65685246838096866817];
LW3_2 = [-0.027531198033021096055 -1.916010297866625578 0.70220025901978844196 8.4270143700945698839 0.00092493896393332371807 0.62910552412053055704 0.56225751760471576191 -40.592501458391460289 -1.0567123736778927245 0.73509462676880621412 -41.760399769139617376 -0.0062966132984655339214 0.12149551864406274437 0.82672236944342014642 -0.0021057698618145975072;-0.00045207873794951070664 -0.084284699665673978863 -0.025655309103261452613 -1.9188542768271430461 -0.0052453187137658082925 -86.600011239166960308 -0.10350401380641299776 -0.056961328138452838044 82.973812309311924196 -0.011845356777081679847 0.32640945648229163778 0.00043205065108231827976 -13.866403329176904435 -0.044665437779255481865 -0.0028811328906409635481;-0.15187017336606115081 -0.03016165248704946733 0.083241653586804451037 0.45728255678961521324 -0.0031274147529882639907 -7.6458504801153495123 0.015780644714278722895 0.082754557708207637789 7.6897693370454849315 0.083087816908937559002 0.12607629545042770003 0.0010661222087953781991 -0.6612999142032136568 0.0069227508827161693405 -0.027401103078615915432;-0.099494336838341274065 2.8375749552748379401 -1.6446406804558613146 -0.95419499139393271214 -1.6626449794984241581 -11.579873386624880638 -0.9259592473610469554 -12.295594199488354192 1.6456519847587942618 -1.5936222379891411549 -11.0759737453655287 0.00091861419377173758518 -0.41954297270509677276 -0.66295826416933745229 -0.0095763535894797646025;-0.41032060187155278141 2.7747236026094772576 -33.132291067703377507 1.9936064069468826254 0.32260976514851463781 -6.2045496443872520942 0.41579242843257047113 -3.1782116823280301965 6.7696246181704768929 -32.455265783443770999 -3.2397152739010843803 -0.067333931805344401922 0.43813144755027955135 1.066144737015576327 -17.699405521570945155;0.056342926634056607338 48.464984148119626184 10.59281269468110942 -11.197346523121790796 0.00093491302803813572213 6.732578363630899787 -17.944899825744897726 1.6314142027917863942 -7.3021479699563176524 11.15822301674768191 12.953489072057084996 -0.14059511541768346388 -0.17997672769054109954 -13.398779773092773127 -0.0024302342934783547658;2.230412065520307241 -11.665974253132484151 -0.7963648304983752535 11.30730154991820946 0.23941459796856687237 66.676267172882560885 4.9179750934930233086 -104.77607178659205545 -70.871953605851132352 -0.78828011190911839723 -108.8447734770000892 1.4944493910913059231 1.8283099500557400408 2.135228438819075425 -3.7621805986740652905;-0.38711254096863922936 -0.076880456604748043103 0.21217800184440394062 1.165600734869802535 -0.007971745142194251732 -19.489049437950900057 0.04022388096651280387 0.21093226725425590073 19.600997344719431226 0.2117859484143303217 0.3213584705584528467 0.0027175026580313768956 -1.6856340160209497547 0.017645958906967140062 -0.06984461024467893131];

% ===== SIMULATION ========

% Dimensions
Q = size(x_sqrt,2); % samples

% Input 1
% no processing

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*x_sqrt);

% Layer 2
a2 = tansig_apply(repmat(b2,1,Q) + LW2_1*a1);

% Layer 3
y_sq = (poslin_apply(repmat(b3,1,Q) + LW3_2*a2)').^2;

% POST-PROCESSING
%[lcts,o2,meth,ac,gal]
y = [-y_sq(:,1:3),-y_sq(:,4)+y_sq(:,5), -y_sq(:,6)+y_sq(:,7), y_sq(:,8)];%[lcts,o2, meth, ac, gal, bio]

% if cell consumes more nutrients than are available: 
ind1 = y(:,1:5) < x_cut(:,1:5);
y(ind1) = x_cut(ind1); 


y1=(y(:,1)<0 | y(:,4)<0 | y(:,5)<0);
%Biomass (flux_bio) is produced only if one of the carbon source is consumed (lactose, acetate or galactose)
ybio=y(:,6).*y1.*(y(:,3)<0);
flux_bio = ybio+(y1==0).*(-vshrinkage);
y=y.*y1;
flux=[y(:,1:5),zeros(Ncells,num_metInert)]; %[lactose, o2, meth, acetate, galactose, bio]

end

% ===== MODULE FUNCTIONS ========

% Linear Positive Transfer Function
function a = poslin_apply(n,~)
  a = max(0,n);
  a(isnan(n)) = nan;
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end
