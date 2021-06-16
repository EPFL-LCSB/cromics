function [flux,flux_bio] = NN15x15_Senterica_pTFA_30K_2rac(vshrinkage,x,Ncells,num_metInert,~,~,rN)

%x=[lcts,o2,meth,ac,gal]
%input for NN: x_cut=[o2,ac,gal]
%Value below which an uptake flux [mmol gDW^-1 h^-1] is negligible because
%the biomass predicted by TFA is zero.
% negl_o=0;   %Oxygen
% negl_a=10^-3;    %Acetate
% negl_g=10^-3;    %Galactose
x_cut = [x(:,2),x(:,4:5).*(x(:,4:5)<-10^-3)]; 

%S. enterica can contain subpopulations. 
%If rN (subpopulation identifier) is not provided then S. enterica is considered as methionine producer mutant 
if nargin ==6
    rN=ones(Ncells,1);
end
x_sqrt = sqrt([-x_cut,rN]');


% ===== NEURAL NETWORK CONSTANTS =====

% Layer 1
b1 = [0.43456359167655295161;0.38010928469245147676;-0.32677956637860666733;0.83975107080931599146;-0.26707379383845225584;-0.56073575957923882385;1.7281367063771375658;-0.70961882701835410359;-0.12082533747361172416;-1.7118273107490364904;-0.2876755821721818096;1.3422657659715240275;0.39177883897439136351;-0.3960539392565950445;-0.47211643232042432849];
IW1_1 = [2.1423139998791911331 -0.0069370188602658755528 -1.5656673699216405726 0.19337357106023725373;5.0774132813054686508 0.02385781631554425089 -7.3921631251033570464 0.20484759802162830211;0.11205145181897718798 0.00074675255389937433316 0.14213436951788221974 1.1604712244254977715;-0.51309690824901577066 0.76782084238256631537 0.68377653004158267525 -0.014040360955373266602;0.00021176869029348465464 9.0573139193120586572e-06 0.056785772639953976948 0.00013916031950447191367;0.12339209282666190048 -3.7179583520687001425e-05 -0.056419968955937239985 0.020055717272219025782;-0.016065284230023853251 0.020874707511331969112 0.86971010397812509307 -0.0037491781595931990452;-4.6498037450293097095 -0.018705889914421541809 6.7638020596241990035 -0.17966835046232623729;-2.0991265196428479101 1.2294021178814231732 2.6054404483090323197 -0.062203905380384942603;-0.73431768594851964949 0.00095694853966604199753 -0.023818084228082923859 -0.024372576149752651214;-1.3761923055141451755 1.4488178067461090492 0.55331871652117492477 -0.0398579734255821716;5.6196935345308585141 -0.0071851661544518642677 -4.6065878592533451297 0.57467662940247732362;2.3425135526348754844 0.0023969090444114571695 -2.5968723424680466572 0.15740186829598201701;-2.1033044911229223928 1.5065465820031147537 1.6662966962511136515 -0.072264946661966922292;-0.735694643134145565 0.014203408486758354809 1.6471470391091767027 -0.089074491573756367502];

% Layer 2
b2 = [17.116646861602561103;-4.5988268061633696959;-4.3131115863652658149;3.4906206449930059499;4.1976404856980886748;4.0474070162943442952;18.829463881944480619;6.3617271178827561329;-2.1912927750274704408;2.9778530855111617903;-30.06694229498425841;-21.984761317094687172;13.858682195340183085;-0.33004835966170548156;42.582421687063309435];
LW2_1 = [0.83045112727704550881 0.73807512970720268175 0.040919913793603004959 0.8007083898335004557 -3.4304305061070681226 -0.85010752982933635025 -13.119782321847745621 6.1059853022329386718 0.029651945080736590021 5.6716627107673946995 0.047150370462841831309 2.1951869836610118369 1.3687799098514368268 -0.075415481198201614332 -0.3697313020249433646;-0.046272665367953859372 0.30321847813348801193 0.051279280137403883411 -0.84077745177343909955 3.274821013779457779 -3.3795830720625277266 16.893853064497339034 -0.019493402592774846521 -0.21714168600671435172 13.142358473287728771 -0.007096827612676123373 -0.035540670826325809761 -0.0091059293893846551815 0.081789096865358609789 0.016940962848857276518;-0.081237614709050651185 0.44308496787474954681 -0.033130418707574008086 -0.0065987276537668075951 -0.18605590236691524342 1.304207862841978649 2.1668922176113563083 0.14565820055459344218 -0.096456856585630926237 -2.6907143902367001154 -0.059595480859390793926 0.094760089623962359684 -0.15353898190514989053 0.019237933374844119749 0.0029846846059883936145;0.0098582976168471378453 -0.062535696907821833324 -0.065443448310976343185 0.57691899087631570708 -2.7196735887221219485 4.5505908278829378943 -11.539801929474204911 0.06381465060687055435 0.10037703784929045969 -10.093029019329888385 -0.0027391513531719795076 0.02511510237208781568 -0.066620992692383598777 -0.050677940438594831774 -0.0941517649174504756;0.013271340413971396618 -0.48771522545746415211 -0.0096639505627861745779 0.7005863421736199248 -1.7607283110072893262 0.64340512658490800479 -7.5441551636607178111 -0.15486286936326390529 0.097915151340008110825 -2.6339155787690762267 0.036937208167045883223 -0.046728772057822547115 0.011337236649884738574 -0.043014981148032621627 -0.11944789327497278009;0.016576461052593850082 -0.45870019526801786292 0.0047151750829211381069 0.60954956068641108491 -0.38520573402273500241 -0.20459900352650600697 -6.6726902493181716736 -0.14265403998294876042 0.091200273974129886234 -1.9705641936205329223 0.039658985313739666956 -0.050951213317449686735 0.02756071745088187197 -0.040113957934512438352 -0.10256524427920266085;0.28620035470182741966 -0.43252696353257291673 0.046622189972995227436 1.0113629345613748267 2.2242153749715543043 -1.963777750013760981 -14.124703885312028362 2.1849416480693499842 1.0071860775754244433 2.1729019593981360181 1.7543946095841558286 -0.94570375169219467182 0.90100905865351110791 -1.168780744614908107 0.3164351476042027711;-0.21597749396194104721 0.061682759666060460457 11.295170868445103096 -0.73067536207932415326 -44.184936643162963321 -13.809774034708048873 -3.8452586278176639212 0.080072662606537614582 -0.061334632394222168839 16.085112046909031847 -0.034190513481709455812 0.056851192092072892759 -0.077153351197519176585 0.087588706314080000404 -0.10750119093822492478;0.30223751057434145029 0.14665830330277992988 0.028751912772469353369 -1.8070966522784717156 3.8991270833473317126 -0.51731649118576983337 8.6005348085325028507 0.083842648891416737822 0.056652150025497963193 4.8984121920622261115 0.049984044798993466918 -0.10442317498328775038 0.42973268776206335895 0.059251645654480031289 0.27642001175070401997;-0.81596164910792212677 1.0103708555608050812 0.0047209980738629605743 0.73121993791300621268 -3.1801613335586673692 0.073228714103068276442 -9.7418827493150406127 2.6213750078559705514 0.15315122817196907823 -5.8723436037568133017 0.118174676365835854 0.56139920217077066145 1.1153765118824421343 -0.16015092534580704919 -0.074078794810052223108;-2.5218689230361062137 1.8928415858010052553 -0.20748385684486672642 -0.97057005375143290404 -2.628088780234085764 8.5593283848246510814 37.756733247926412389 -3.156912515945869746 -3.4547445133057963496 4.4047153831215890563 -6.0066290089965495014 2.1331828897663296374 -1.5509405552902237879 6.697194321653367588 -1.8056014129967787873;-3.8844344776316619416 5.4134811298105320887 0.042310866262060466136 -0.65926200389039690819 8.5667987538368723932 -2.4811152030981231675 4.5414475426888447629 4.1530023083160667596 -0.026259869307900012697 -16.275104724727963657 -0.028762261324929816464 5.3828457492506158744 0.67750174881819003048 0.054069571749004959049 0.30821240439493435836;0.79268668522689000877 -0.62973771439490244806 0.029715298330674901756 0.91076143081699334481 3.7917646231645192323 -1.4216291351475296434 -9.2725500841699144416 0.45864861066610790319 0.66395241843431740403 2.6968058759118886947 1.0967152407239675327 -1.5281281276519804724 -0.0071034257627956920791 -1.165315676996252936 0.050419093289587657569;-7.0014181838335185328e-05 9.5100687139683973094e-05 0.00012568522669586431055 -3.596159915166415147e-05 -0.60577892890350470978 -0.0027935150809036009799 -0.0034708174916997527014 3.2841243125373823137e-05 -2.5052411890139640778e-05 0.0035209641576640759303 -5.2501812492222913007e-06 -7.7172507290522851125e-07 2.9147135103356379757e-05 2.8323319520428663068e-05 -4.5461893652035276702e-05;1.7885052754303276679 0.87206399141612644943 -0.029372199371165322312 1.8912498992305299961 3.7354403631909374184 0.57448711301151667019 -11.59638168539146541 1.0168059763817227648 -0.41664625569314961417 0.31521418903694348179 -0.56127968919184323404 -28.817729034740739991 -4.4192204933778045373 0.60452335117642885809 -0.49094291659058869604];

% Layer 3
b3 = [-6.5012354449841236459;-2.1862647537163160116;2.7983430730744456305;-5.6682402542096426146;-5.4972742073859688361;-1.2674566488058456315];
LW3_2 = [0.014216192544294408914 2.565946988861262934 0.67865124979955948614 2.592288261526468407 9.9857866463748603536 -9.8666630636035836943 -0.94858718510981199223 -0.051185988767906101526 0.62859389043604096958 0.025245086834202141307 0.091109092108481615213 0.0062750811236395673959 2.5861804183766792242 -33.838469506221422023 -0.0040350159073331073878;0.008519467110550765776 0.26419399632281498347 0.20496503635728766102 0.15215243492638103051 0.94621327453987447154 -0.7851793961536368549 -0.092847490755304407095 1.276928205315854381 0.1015807095805478788 0.0035370473633303946108 0.0087240875575999091152 0.0067105156571922851966 0.2547621339889784009 -5.473664234493599956 0.0022245446231041496148;2.2642195299971898947 0.10200649850726892232 2.4423313815518903347 1.2805839946218808834 19.436887982692422128 -9.2534243925649057871 -1.1776737780756620921 -0.024733811204446471549 -1.1724816436139666731 -10.959326457626350049 0.11697180215407522452 1.9378266876893472936 3.4250184463892869857 24.619986419542744471 -1.8480954490797740775;2.4090396608985549243 11.643660174979373778 9.5388683462983170358 13.558403225000514425 -24.649361794213611176 37.155799252775786101 4.6172200794150954906 0.035256382772809173598 1.6428347072799693862 0.59600563799556927069 0.15471437074035845893 1.2491847942137850325 0.54804255468938922746 16.267944368074683581 0.1374294487249406449;-0.0032158489167102878581 0.014171407859912029184 -0.020170556116373857675 0.01153791260618410898 -0.10991819005498947026 0.099325705257839166928 -0.0071395731009420253244 -0.00052760825711530717503 -0.010665401032316836408 -0.00045229194555595131008 0.00054323325049143321355 -0.0011772538133633073283 0.011224134536999136325 -31.270087772638188284 -0.00079566490328302852484;0.011563121101539258762 0.37648767942216260218 0.28915186485855731702 0.21758444192406992124 1.3407376222425688717 -1.1147494896528353792 -0.13051268186610565092 -0.018648614401018236836 0.14401023216695305051 0.0053135938654130269662 0.012402655056235336109 0.0089044434104908511129 0.36033771445067575012 -7.7398130051722873191 0.0034675222613220514634];

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
y = [-y_sq(:,1), -y_sq(:,3)+y_sq(:,4), -y_sq(:,5), y_sq(:,2), y_sq(:,6)];%[o2, ac, gal, meth, bio]

% if cell consumes more nutrients than are available: 
ind1 = y(:,1:3) < x_cut(:,1:3);
y(ind1) = x_cut(ind1); %[o2, ac, gal]


y1=(y(:,2)<0 | y(:,3)<0);
%Biomass (flux_bio) is produced only if one of the carbon source is consumed (acetate or galactose)
flux_bio = y(:,5).*y1+(y1==0).*(-vshrinkage);
y=y.*y1;
flux=[zeros(Ncells,1),y(:,1),y(:,4),y(:,2:3),zeros(Ncells,num_metInert)]; %[lactose, o2, meth, acetate, galactose, bio]

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