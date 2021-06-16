function [flux,flux_bio] = NN15x15_Ecoli_glc_o2_30K_pTFA(vshrinkage,x,Ncells,num_metInert,~,~)


%Value below which an uptake flux is negligible
negl_bio=6e-4;   %Biomass [h^-1]. Assumed

%input=[glucose, oxygen]
x_cut = [(x(:,1)<-0.14).*x(:,1), x(:,2) ]; 
x_sqrt = sqrt(-x_cut');

% ===== NEURAL NETWORK CONSTANTS =====

% Layer 1
b1 = [0.61751454008084216696;-60.140582435745308487;0.75026208094269619675;2.062429820930004265;0.40898529991284848784;0.49910610518010106729;-0.1916094063076893006;1.6770586267550235959;-13.182991405797782392;-1.4710884776405657171;-32.638505619111384704;3.8071501546729029641;-0.68680582682260482574;60.121155687893136133;1.4129966293774613995];
IW1_1 = [0.20245380027974732573 -7.2548586898892386543;51.233493383768525575 44.949097460247543268;-0.60376812880690122753 -0.0037370135615584277634;-3.9128015020668671653 0.011445416278353633469;0.60052286295031820362 0.0008550397129562729414;-1.0114400422390834589 1.1007798705365072767;2.5187990637421839857 -2.2671100353769046976;11.735866800525275977 -15.727799716158951782;37.491675236740441335 -0.021692426255494168252;-7.2518952804214622532 6.0317806557213504703;29.240271295525026574 2.9945915195126855401;-11.113118354264887699 -0.016541732488440933496;0.47105773282001023317 -0.9682456537755683712;-51.221780217930188428 -45.870141813962611366;-0.38789541217194001277 0.00087587971721624590087];

% Layer 2
b2 = [0.36752548178655791711;-17.619075575776577125;-4.714875352230557759;-7.3823722436776950673;-8.3228017906440570073;2.2605991840360952772;0.29210272637958428765;-19.087941663910363843;6.4784456338582696588;8.8265189662276508642;0.081454499372989805361;-12.729882764031991371;-1.7780721695258514536;-0.88969246149691305625;-14.307237651311469406];
LW2_1 = [-0.0072171404365577130832 1.5323081369692788112 0.0022381636172191417267 -0.0320677127804149531 0.64144389460966999561 0.019432754761548692163 -0.082855443294848077351 0.026525701744530623405 -0.08033974388709448633 0.13627336542550300136 -0.0021867851134615978737 0.31320512482651124975 -0.033269228231729008893 1.663019454713806411 -0.20694800364102489132;13.770400670463347481 31.670131738472679928 -2.3098217844695310141 -0.098624109656272634927 -4.7853007842298067587 -9.5031352648295364105 -4.9194965884789123578 12.022298238271758919 -5.1025760594523692504 -4.3842549641451551423 -0.16033302985980615873 -2.6374873929670461514 18.740667728434868167 -43.424299008973413549 -3.2414729438053098143;0.37726746607646843845 5.2403188667433555636 -5.007359857835867345 -7.162933388580819738 0.85137486708861165319 4.7158425190600778976 -3.6309245785499912884 -3.6792551097184480469 15.161534084112940235 3.179000528621287458 21.892076318068905749 -7.5450695059592591818 -4.4765495968182316133 -5.746582448094805784 -5.0223830853088085036;2.4004462806359980753 -4.1648147261014072384 2.2069278102612370773 0.045365699363804967703 4.5961255737593891268 7.6623216753250140698 4.6303869091799487734 -1.1251570723969632759 -3.3798505804778700146 4.3803607591401911847 0.13606226668215140596 4.5133777954496974161 -15.727355226686361434 4.7500337104774636998 3.2299214622856378121;1.460071905199535669 -6.5524563184639870173 1.7033498166941070107 0.042767242536080087023 0.92387315024114857742 2.6478778589197480287 1.8138943195689163357 2.3944749440747115088 -7.387154728788440039 3.9679127453816911419 0.15357427215297170253 2.9663331620780688525 -28.385497235731186549 2.6554866088173350036 2.8973005118188677187;-0.029108685555372874887 -3.717521126333843462 0.16602436709994100594 3.8933801750592249569 0.59318850173527426506 0.40323562750980235236 -0.490021470087717681 -0.0089325381494076549344 1.1330865247554160469 0.1379953988435914658 0.0029503170946727310699 -1.4809604066768631281 -0.34306871567265695422 -4.0690076541261808174 0.27452342483389297056;-0.0063932878552647245951 -0.2639492751721969066 0.016199118641676103075 -0.016425780071878698374 0.58815365906363126225 0.067377519692560478015 -0.11052583545056167413 -0.25879926583650025851 -0.10357611102585560781 0.17757411482093374278 -0.0020478760725517013153 0.61792778595534747499 -0.062598422638507969507 -0.20011337638007542505 -0.12722764533294042555;39.877484912440145592 43.621711591870294455 1.4079079468862178981 -0.11225684669979751895 -5.2985809814916295934 -13.927332048996287739 1.214933811069260905 0.024259783853194046704 -18.581344129891952832 -5.8081900141646203295 0.00066797548387624995463 11.837634781397190054 0.44170714133017729974 -44.491242789191474571 6.1183326987797270746;2.071365318098361552 0.42991184652784869558 -0.53893589620589410494 0.25376195685828495385 -4.3480164230476061604 4.2460834231412567519 2.9570906615413141694 1.5224178492144833097 1.9831748584609774078 -0.3219779952347950891 0.031291896800777506404 3.220697040682388046 3.2024322879473827719 -2.987015458833639503 0.52944797783954611337;4.1715985240048292226 -13.218026272308154745 15.76930416424888115 20.518717060938477914 -3.5206471850777494303 -5.4151107169103118721 5.7640037335037392907 2.4699401042270618944 -24.428975589495749432 -16.823465290285344764 -46.998129568777009979 13.319357331897750285 14.639827242608546953 6.7271220872271584312 10.307650246284309148;-0.00025110940200897050294 0.036154069007174292893 -0.028287925123802300231 -0.00070305725946008328867 0.07009707493661226918 -0.00061936003778801848856 -0.002064142122229020021 -0.02425932854752602813 -0.0056151210646938750845 0.0035353907693460176823 -4.6873336772485227678e-05 0.007133035981476220852 3.0872624221344297207e-05 0.040275159417759227864 -0.11626281384551685749;16.050364920927904677 14.462595432406395446 0.10287353300997852112 -0.075993639004765037304 0.81306068255737606876 0.34036230361181424309 -0.16783426634710527092 1.3239521832078857688 -1.2458772931926322958 0.27548377461902490149 -0.0027659168525142369295 -0.37957553370252000713 0.89859609365757442401 -16.403357423451751629 0.096264872124410913456;0.0028673780091982236339 2.5396355752243082726 -0.22885195072105335434 0.14552652861112863603 -0.95871518080193507405 -0.2869877360445011738 0.37499425989752266375 -0.08107928373003112088 0.66088267677725009985 -0.33914850599254864338 0.0055861346477929767901 -2.2297838587168921087 0.38610944784083167658 2.7141753210704746557 -0.26674462353372135892;0.017465233077386066451 -0.71218472192721105074 -0.36381343398119631027 4.1338968204264761042 -0.17914498024674774257 -0.50974593007835333758 2.1547050901218534058 0.014846179987354305232 -1.8905819861201009768 -0.49618892230532618548 -0.016460575035632064866 -3.0073532467416845826 0.24814885621506829749 -0.36048017825161432492 -0.34513559170318414138;16.652446008139275335 14.941442911210749855 0.091043380540224566611 -0.081532583475807107742 0.32102443528125712557 0.091636619161843804737 -0.096245067155867267128 1.2724470149002542474 -0.90903796779862289235 0.064489782081330757069 -0.0020744014410828482201 -1.0997922710772716215 0.74071409285043821047 -17.057570943133448793 0.18785087985330659044];

% Layer 3
b3 = [-5.4759203144580856915;-12.327105939299373816;-0.37052937798055513019;-2.4199833068602174535];
LW3_2 = [-0.79671092798627807241 -0.6273833906594268317 -2.020733725464588737 0.096128126862394552843 -0.71812991720339569568 -0.040738230693549218331 0.38606313369593359131 0.63840431452218593389 -0.0029558065109108570823 -8.1729486563430331358 20.993361503706744031 -0.21451066078890951294 -0.022930535873443930184 -0.091111548069105693926 0.20294688806526997582;-5.1782009426126505147 0.89498182454687236209 -5.7555718894787588624 -0.18849165540359799631 1.1040412460032138497 0.84944448481797707462 1.41584806397979035 1.1174389199714087617 -0.13707458325100885244 -11.530929947732312613 30.880801258829023936 0.63276722316506983645 -1.1078218294749193618 -6.7010759965169128805 -0.60805806686371866032;-13.327299083237589272 -7.6346168516435710671 -3.5076401957044089919 1.4536001534766447296 -9.1689011498354808793 -1.4553263521089510046 25.56247480391609983 6.3512120215709577664 0.55385727041565113726 -12.654331227228251322 14.729035789423379654 -11.387236427125058569 2.0848186280037550588 0.25059132033241299231 12.667002266104262276;0.20693570226934956957 -0.021361894065905561585 0.82368083855808282845 -0.0029703340106094970434 -0.015754351406341723929 -1.1155250886085181783 -0.44579782090758607316 0.033690694448152674889 -0.0025095829136219144524 -0.22867428967198266676 6.6740639681641109604 -0.06303126844761636205 -0.291268124685951868 -2.3498453315992784418 0.050280084499835009382];

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
y_sqrt = (poslin_apply(repmat(b3,1,Q) + LW3_2*a2)').^2;


% POST-PROCESSING

y = [-y_sqrt(:,1:2), y_sqrt(:,3) , y_sqrt(:,4).*(y_sqrt(:,1)>0)];

% if cell consumes more nutrients than are available: 
ind1 = y(:,1:2) < x_cut(:,1:2);
y(ind1) = x_cut(ind1);

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
