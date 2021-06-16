function [flux,flux_bio] = NN15x15_Ecoli_glc_o2_ac_30K_paTFA(vshrinkage,x,Ncells,num_metInert,~,~)

%Value below which an uptake flux is negligible
negl_bio=6e-4;   %Biomass [h^-1]. Assumed

x=-x;
aa=(x(:,3)+3.75*x(:,1))>=0.6; %Based on TFA evaluations, 0.6 is the minimum amount of glucose and/or acetate required to produce a biomass flux greater than zero.
x_cut = [aa.*x(:,1), x(:,2) , aa.*x(:,3)]; %[glc,o2,ac]
x1 = sqrt(x_cut');

% ===== NEURAL NETWORK CONSTANTS =====

% Layer 1
b1 = [-1.8935743037428249824;-21.9695772146962085;26.226015660345975533;-0.43158860615948635431;1.3169294835771347962;1.7697848357666281505;0.51118973145922297352;0.86005863042316987865;0.19714187964288026889;-2.3833672458794170623;1.1363943298835181839;-0.22480052040231404686;-2.8042364347988502082;-1.1918221861636191239;0.56643935586224625212];
IW1_1 = [-2.3360893064905186023 2.3002115413000638888 -1.5115051271591075643;3.954295572631525868 3.0777263476598522018 -0.0026681186982883231366;-2.6241041310439858769 -5.7449689387806417429 0.0030778328034045082488;0.34415913727889557716 -0.22665422519184050265 0.00022676806669735135821;0.44464675693549554625 -0.36743896199816877823 3.4439488728882044328e-05;1.4825873389202319785 -1.5289040052860780605 1.2446660406371772289;0.52372476372173659698 -0.57590685748179903936 0.00024411401494541223184;-0.36134141200702657448 0.09551032661119969458 -7.5913008137643272349e-05;3.7219343907053317189 -3.6840774885108524828 0.0013033252338981895883;-23.589537428604586466 18.463321167743632856 0.00050240999159969872145;3.9542112000998779209 -3.8406132932975589611 2.3216787152096189928;0.53523709743258096605 -1.053455003226742015 0.0012867026175319448787;-1.2884376549005587265 1.3288466568111620969 -1.0927025304772719494;-0.43964723969683638938 0.32334719656132560051 -0.00014808184144144394072;1.7990614749943829587 -0.31462013813196115786 0.0039016897581533650265];

% Layer 2
b2 = [-0.39678194366376179358;0.5458275021497344115;2.3343331256115331662;-5.498586286918716759;-8.9988901863388619518;-0.42566397338174982723;0.30122224638780514816;0.57884469059749554809;-3.7697467407191833644;64.296798840227708638;1.4057836030610673106;-137.58810985681913053;-0.17991775861132697223;-87.823878077628236838;88.318434793923415782];
LW2_1 = [0.0076966053611942615922 -0.59928551472822766399 0.010423588825242297365 0.16127118914349353185 5.9644532082054571021 0.041291452473268590306 -0.094546367436755485159 0.639767307471933111 0.048031805160444772385 0.028795900473423514943 -0.00096762972118756909734 0.038389665855935353744 0.060192299164926500732 6.9424204746379585274 -0.066255508529067488888;3.1325509288433770472 2.2975506719085596607 -0.010856407346147596907 -1.9898283150804829678 8.8749929470350465976 3.4642126398618073857 -0.72342444886663470616 -0.13276840595906108655 -1.8496397157952402246 -0.3121047401795723486 -0.18128789032182623941 -0.233776702382713919 0.16459919679217097088 8.5926404338619359891 -0.052987434107237191272;-0.45628809567987466123 -0.2028984995852662343 -0.076405854877075463083 33.746044758668503505 -37.772534447357948295 0.2344415822357810153 -7.9848961335659955907 20.896862421839255575 -0.10125289988102938177 -0.045026559039816227425 -0.049639497284977872338 -4.351115035986066637 -0.54458269706850637171 -39.157313331856293814 -0.12907124398552027866;-0.19481484667020895762 0.15357171304414307667 -0.13069164629770654851 3.9679552605125953946 6.535556016171445215 -0.22561487801061178016 -2.8116139776706976505 2.3309180024009363841 0.27958328621821210147 -0.0073807691669643456928 -0.03030475733332736038 0.90047734700845882827 -1.6477256795258599098 5.3715465423633519393 1.8983951474909168677;-0.8175002083742090031 1.352859571199519273 0.21785773452040060882 1.9773223383276301579 -9.7003620706870581358 -2.1556345837185761916 0.32125836611259844799 1.3456196593085678614 -1.3474790381643060311 -0.34496524104708586878 0.3713218210997494273 0.32938975361980182832 1.512921453149952411 -24.226929719251170781 -1.2265888464126331048;0.14530101255096805613 0.5973757069050225077 -0.019162247574856291349 1.4560895795908779249 -8.3552941029203218193 -0.096382130657435519105 -0.19214021083056867512 0.14890862049811132151 -0.011359950181415350359 -0.03522054665254671646 0.040009737472843316664 0.087047334798832184943 0.36575601157559939347 -9.8553802557873488865 0.79174148906577745954;4.8816053458205104221 1.7316801893008795243 0.17718893165381308097 3.5465405969786085727 -6.5159400886637879324 1.4630073408409700608 -1.1878857046855473278 3.0004411976446307975 -1.0230781510324344907 -0.3432737617568561439 -0.11985000930316801193 -0.26763922273612961167 -3.6623411217652277472 -8.2946132523132050096 -0.17522159538019271841;0.025818942063842983958 -0.75960333311722638694 0.013582906237010840267 -0.0074397718347921715579 8.8442951509983398495 0.06716976807805861982 -0.2527006026599791233 -1.0083018673945149946 -0.0023322396222364791957 0.059443249259332420209 -0.014362708985323434027 -0.047819146766597361853 0.14330351995982526514 10.406636015880200929 -0.074596207216477250501;0.016668427775226857951 -0.53365786198202835866 0.027606509906443543051 -2.9580142915101137469 3.7846036549848518504 0.13863539822632481302 0.097325044715112968285 -2.1020562318573685623 0.037598102226217677835 -0.03897356222703637213 0.013652237799764283671 -2.2286720999877434046 -0.034499135150119625282 4.653346767483309776 0.043897733916474994398;-0.035469124753165305608 -5.5050391385115169385 2.3102660741503058262 -4.9738698867455495289 -126.35055980587294755 -7.990449070140392962 5.7965664325811374624 -11.911099476878655778 3.2639737005948772364 1.1044678960139044932 0.45476512302703780977 -1.1925061755820995923 -40.774509757134062227 -4.9063297297854111534 4.0779847817210441718;0.12403870907881810182 0.39361033428737995932 -0.049814255636348017675 7.4224837489226542786 0.73873249625663384155 -0.20161762765786467311 -1.1965198586331970976 3.9263049646122749614 -0.15670978695134346603 0.0046127913359565592133 0.044536243182492982196 2.2306545865432370057 0.74606676016034323418 -4.4205703620513325802 -0.79606956836818920653;-6.3221627407944831489 10.401408429452446924 -1.5000815549096551127 -3.3518569056108988846 -203.39126016666759256 -28.647565640215248095 1.6145243744843502665 13.754486526800578261 0.27691676830661487996 1.3374956950545433543 1.0779348883841068929 -6.6290707236463513041 -131.39709908313957953 -245.88611403327465155 -8.5451605093182880779;0.53218437603879220621 1.6970976062916844906 0.041242755069744604923 3.1696171732726905468 -22.273973889385739966 0.029908809539988399223 0.13439521838486642724 2.6170724572263281971 -1.2162391063885829023 -0.33116541459496751099 -0.2404580843362374154 0.218951838444423591 0.5571513133986725208 -26.112247886464253099 0.087402864417179509937;-5.1649254738357486261 10.273506276241036517 -0.5883511475171313565 6.4578094068948255213 9.5175038271127512246 -26.704092531070564576 -0.60016137701917160463 3.6321040789372260882 -0.32864965667682843886 0.58551607155035434449 1.1223693446024693987 0.25353516194896857927 -128.18835113472479748 18.685919685467048623 -2.0353095561884315678;5.189963596801694834 -10.314782281228497496 0.57742261669631467225 -6.2212257026686987871 -8.5320733478102539493 26.759517588846701841 0.56206279260623315253 -3.6387268195689643768 0.32067124764982773977 -0.58509686973663332932 -1.1193080906044567957 -0.13853307364233291565 128.43664959122099845 -17.463982213875784311 2.0896311938967868649];

% Layer 3
b3 = [2.796349007397477493;28.700610072172022313;-44.993941225677325235;74.362173209488318548;-49.319868629658643044];
LW3_2 = [-4.8238823768748186893 -0.071888696798133866483 -2.6036647099909466796 -0.92881077494300534614 0.10360102235502922918 -0.098394253831061276316 -0.13017316110106741389 2.3766341675704940606 0.41671457123916233467 0.0034580704217355611244 -0.053653027578795177421 -0.015613402573730309972 -0.060496569296800833948 0.07210931504693024463 0.074980301691402073683;-4.9368080625764472202 0.12003190629348851892 -21.573668060223116782 1.0621731301426631244 -0.19955926364491904934 -0.059074527585686681386 0.94176268650146077732 1.249819717823718479 6.8434100651883547073 0.024733871555722576618 0.33843570423600294328 0.0033860844146036362876 -0.88911692002242881383 -6.8622755577819178896 -6.8537102899496167296;-0.62644030287774121746 -5.6613186149802325886 -50.325866478525270509 8.4214080739679673115 -11.744341912122475691 -1.0620748658466319103 2.1762112060347376108 -2.8874133254446641139 8.276203401110967306 102.27583988060096942 -0.29654623963069137726 0.039321240877837039529 -0.728324053298657903 -8.4237324319356705615 -8.4073966532818271702;-18.414343362402508575 -0.4115547742942020637 -24.632966253391710865 55.022002786336841496 0.22527556982733035196 -8.5156321185490462256 -2.2793408805906483217 1.117489660063070378 7.9156096761846388432 0.12933219723258770895 6.8202109050780803656 -0.93202076322269999942 2.1740886449105909506 -50.92007206672534636 -51.652065797977876116;-1.2785730777894284138 0.0088374941540650345934 51.10036853367044074 0.64210060040318517327 -0.018187279805981920733 0.05934550870399756678 0.17029196520329378806 0.54299073450719304024 1.1849544593708927209 0.0060413984670681947914 0.076429483468516093425 -0.0026050128191708633862 -0.16539988718487380615 -0.98821637436875719995 -0.98741924148008974882];

% ===== SIMULATION ========

% Dimensions
Q = size(x1,2); % samples

% Input 1
% no processing

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*x1);

% Layer 2
a2 = tansig_apply(repmat(b2,1,Q) + LW2_1*a1);

% Layer 3
y_sq = (poslin_apply(repmat(b3,1,Q) + LW3_2*a2)').^2;

% POST-PROCESSING
%[glc,o2,ac,bio]
y = [-y_sq(:,1:2), -y_sq(:,3)+y_sq(:,4), y_sq(:,5)];

% if cell consumes more nutrients than are available: 
ind1 = y(:,1:3) < -x_cut(:,1:3);
y(ind1) = -x_cut(ind1);

y(:,4)=y(:,4).*(y(:,1)<0 | y(:,3)<0);
flux_bio=(y(:,4)<negl_bio).*(-vshrinkage)+(y(:,4)>=negl_bio).*y(:,4);
flux=[y(:,1:3).*(flux_bio>0),zeros(Ncells,num_metInert)];

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