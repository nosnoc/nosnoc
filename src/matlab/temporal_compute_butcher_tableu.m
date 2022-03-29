%%  Legendre collocation points
   legendre_points1  = [ 0.50000000000000000000 ];
   legendre_points2  = [ 0.21132486540518713447, 0.78867513459481286553 ];
   legendre_points3  = [ 0.11270166537925824235, 0.50000000000000000000, 0.88729833462074170214];
   legendre_points4  = [ 0.06943184420297354720, 0.33000947820757187134, 0.66999052179242823968, 0.93056815579702623076 ];
   legendre_points5  = [0.04691007703066807366, 0.23076534494715861268, 0.49999999999999994449, 0.76923465505284149835, 0.95308992296933192634 ];
   legendre_points6  = [ 0.03376524289842347537, 0.16939530676686742616, 0.38069040695840172805, 0.61930959304159849399, 0.83060469323313235179, 0.96623475710157580298 ];
   legendre_points7  = [ 0.02544604382862047931, 0.12923440720030288098, 0.29707742431130129690, 0.50000000000000000000, 0.70292257568869853657, 0.87076559279969734106, 0.97455395617137896558 ];
   legendre_points8  =  [ 0.01985507175123157886, 0.10166676129318691357, 0.23723379504183561561, 0.40828267875217505445, 0.59171732124782483453, 0.76276620495816449541, 0.89833323870681347501, 0.98014492824876797705 ];
   legendre_points9  =  [ 0.01591988024618706810, 0.08198444633668211523,   0.19331428364970504319, 0.33787328829809543107, 0.49999999999999988898,  0.66212671170190451342, 0.80668571635029517886, 0.91801555366331766272, 0.98408011975381259884 ];
    legendre_points =  {legendre_points1, legendre_points2, legendre_points3, legendre_points4, legendre_points5, legendre_points6, legendre_points7, legendre_points8, legendre_points9};
%%  Radau collocation points
   radau_points1 = [1.00000000000000000000];
   radau_points2 = [0.33333333333333337034, 1.00000000000000000000];
   radau_points3 = [ 0.15505102572168222297, 0.64494897427831787695,1.00000000000000000000];
   radau_points4 = [0.08858795951270420632, 0.40946686444073465694,0.78765946176084700170, 1.00000000000000000000];
   radau_points5 = [0.05710419611451822419, 0.27684301363812369168,0.58359043236891683382, 0.86024013565621926247, 1.00000000000000000000]
   radau_points6 = [0.03980985705146905529, 0.19801341787360787761, 0.43797481024738621480, 0.69546427335363603106, 0.90146491420117347282,1.00000000000000000000];
   radau_points7 = [0.02931642715978521885, 0.14807859966848435640, 0.33698469028115418666, 0.55867151877155019069, 0.76923386203005450490,  0.92694567131974103802, 1.00000000000000000000 ];
   radau_points8 = [ 0.02247938643871305597, 0.11467905316090415413, 0.26578982278458951338, 0.45284637366944457959, 0.64737528288683043876, 0.81975930826310761113, 0.94373743946307731001, 1.00000000000000000000];
   radau_points9 = [ 0.01777991514736393386, 0.09132360789979432347, 0.21430847939563035798, 0.37193216458327238438, 0.54518668480342658000, 0.71317524285556954666, 0.85563374295785443735, 0.95536604471003006012, 1.00000000000000000000 ];
  radau_points =  {radau_points1, radau_points2, radau_points3, radau_points4, radau_points5, radau_points6, radau_points7, radau_points8, radau_points9};

%%
tau_col = legendre_points9;
n_s = length(tau_col)
y = sym('y', [n_s 1]);
syms tau

lagrange_poly = 0;
lagrange_poly_integral = 0;
for ii = 1:n_s
    term_ii = y(ii);
    for jj = 1:n_s
        if jj~=ii
            term_ii = term_ii.*((tau-tau_col(jj))./(tau_col(ii)-tau_col(jj)));
        end
    end
     lagrange_poly = lagrange_poly + term_ii;    
end

% lagrange_poly = simplify(lagrange_poly);
lagrange_poly_integral  = int(lagrange_poly,tau);
lagrange_poly_fun = matlabFunction(lagrange_poly);
lagrange_poly_integral_fun  = matlabFunction(lagrange_poly_integral);

%% Butcher tableu
c = tau_col;
A = [];
b = [];
eval_str = [];
for ii = 1:n_s
    eval_str = [eval_str,  'y(' num2str(ii) '),'] ;
end
eval_str(end) = [];
for ii = 1:n_s
    a = [];
    y = zeros(n_s,1);
    y(ii) = 1;
    for jj = 1:n_s
        a  = [a, eval(['lagrange_poly_integral_fun(c(jj), ' eval_str ')' ])];
    end
    A = [A;a];
    b = [b,eval(['lagrange_poly_integral_fun(1, ' eval_str ')' ])];
end
sum(b);
A = A';

% [A_irk,b_irk,c_irk,order] = generatre_butcher_tableu(n_s,'legendre');
% % 
% % 
% error_A = norm(A-A_irk)
% error_b = norm(b-b_irk)
% error_c = norm(c-c_irk)

% A
% b
% c

%% print results
clc
str_c = ['c = [' num2str(c,16) '];\n'];
str_b = ['b = [' num2str(b,16) '];\n'];
str_A = ['A = [' num2str(A(1,:),16) '; \n'];
for ii = 2:n_s 
    str_A =[str_A num2str(A(ii,:),16) '; \n' ];
end
str_A(end-3:end) = [];
str_A = [str_A '];\n'];

% str_A = ['A = [' num2str(A(1,:),16) '; '];
% for ii = 2:n_s 
%     str_A =[str_A num2str(A(ii,:),16) '; ' ];
% end
% str_A(end-3:end) = [];
% str_A = [str_A '];'];
% 
% eval(str_A)
% error_A = norm(A-A_irk)

fprintf(str_A)
fprintf(str_b)
fprintf(str_c)


