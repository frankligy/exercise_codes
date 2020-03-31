using Distributions;
using DataFrames;
using RDatasets;
using Gadfly;
using Cairo;
using Fontconfig;
using CSV;

function E_step(x,mu_a,mu_b,mu_c,sigma_a,sigma_b,sigma_c,pi_a,pi_b,pi_c)
    numerator_a = zeros(size(x,1));
    numerator_b = zeros(size(x,1));
    numerator_c = zeros(size(x,1));
    denominator = zeros(size(x,1));
    post_a = zeros(size(x,1));
    post_b = zeros(size(x,1));
    post_c = zeros(size(x,1));
    for i=1:size(x,1)
        numerator_a[i] = pi_a .*pdf(MvNormal(mu_a,sigma_a),x[i,:]);
        numerator_b[i] = pi_b .*pdf(MvNormal(mu_b,sigma_b),x[i,:]);
        numerator_c[i] = pi_c .*pdf(MvNormal(mu_c,sigma_c),x[i,:]);
        denominator[i] = pi_a .*pdf(MvNormal(mu_a,sigma_a),x[i,:]) + pi_b .*pdf(MvNormal(mu_b,sigma_b),x[i,:]) + pi_c .*pdf(MvNormal(mu_c,sigma_c),x[i,:]);
        post_a[i] = numerator_a[i] ./denominator[i];
        post_b[i] = numerator_b[i] ./denominator[i];
        post_c[i] = numerator_c[i] ./denominator[i];
end
    return post_a,post_b,post_c;
end

function M_step(x,post_a,post_b,post_c)
    mu_a = sum(post_a.*x,1)./sum(post_a);
    mu_a = Vector(mu_a[:]);
    mu_b = sum(post_b.*x,1)./sum(post_b);
    mu_b = Vector(mu_b[:]);
    mu_c = sum(post_c.*x,1)./sum(post_c);
    mu_c = Vector(mu_c[:]);
    sigma_a = round.((post_a.*(x.-mu_a'))'*(x.-mu_a')
              /sum(post_a),5);
    sigma_b = round.((post_b.*(x.-mu_b'))'*(x.-mu_b')
              /sum(post_b),5);    
    sigma_c = round.((post_c.*(x.-mu_c'))'*(x.-mu_c')
              /sum(post_c),5);      
    pi_a = sum(post_a)/size(x,1);
    pi_b = sum(post_b)/size(x,1);
    pi_c = sum(post_c)/size(x,1);
    return mu_a,mu_b,mu_c,sigma_a,sigma_b,sigma_c,pi_a,pi_b,pi_c
end

function EM(x,mu_a,mu_b,mu_c,sigma_a,sigma_b,sigma_c,pi_a,pi_b,pi_c)
    maxIter = 1000;
    for i=1:maxIter
        print(i,"\n");
        post_a,post_b,post_c = E_step(x,mu_a,mu_b,mu_c,sigma_a,sigma_b,sigma_c,pi_a,pi_b,pi_c); print([post_a post_b post_c],"\n");
        mu_a_new, mu_b_new, mu_c_new, sigma_a_new,sigma_b_new,sigma_c_new,pi_a_new,pi_b_new,pi_c_new = M_step(x,post_a,post_b,post_c);
            print(mu_a_new," ",mu_b_new," ",mu_c,"\n");
            print(sigma_a_new," ",sigma_b_new," ",sigma_c_new,"\n");
        if(sum(abs.(mu_a-mu_a_new))<0.001 && sum(abs.(mu_b-mu_b_new))<0.001 && sum(abs.(mu_c-mu_c_new))<0.001
                && sum(abs.(sigma_a-sigma_a_new))<0.001 && sum(abs.(sigma_b-sigma_b_new))< 0.001 && sum(abs.(sigma_c-sigma_c_new))< 0.001)
            break;
        end;
        mu_a = mu_a_new; mu_b = mu_b_new; mu_c = mu_c_new;
        sigma_a = sigma_a_new; sigma_b = sigma_b_new; sigma_c = sigma_c_new
        pi_a = pi_a_new; pi_b = pi_b_new; pi_c = pi_c_new;
end
    return mu_a,mu_b,mu_c,sigma_a,sigma_b,sigma_c,pi_a,pi_b,pi_c
end

iris = dataset("datasets","iris");
x = iris[:,[:PetalLength,:PetalWidth]]
x = convert(Array,x);
mu_a=[1, 0.5];
mu_b=[2, 1];
mu_c=[4, 1.5];
sigma_a = [10.0 0.0; 0.0 10.0];
sigma_b = [10.0 0.0; 0.0 10.0];
sigma_c = [10.0 0.0; 0.0 10.0];
pi_a = 0.4;
pi_b = 0.3;
pi_c = 0.3;
mu_a, mu_b,mu_c, sigma_a, sigma_b,sigma_c,pi_a,pi_b,pi_c = EM(x,mu_a,mu_b,mu_c,sigma_a,sigma_b,sigma_c,pi_a,pi_b,pi_c);

print(mu_a,mu_b,mu_c,"\n")
print(sigma_a,sigma_b,sigma_c,"\n")
print([pi_a pi_b pi_c])

##plot the result
data = dataset("datasets","iris");
data_mat_a = data[find(data[:Species].=="setosa"),[:PetalLength,:PetalWidth]];
data_mat_b = data[find(data[:Species].=="versicolor"),[:PetalLength,:PetalWidth]];
data_mat_c = data[find(data[:Species].=="virginica"),[:PetalLength,:PetalWidth]];
nrows_a = size(data_mat_a,1);
nrows_b = size(data_mat_b,1);
nrows_c = size(data_mat_c,1);
#Estimate these using EM for MV Gaussian approach
mean_vec_a = vec([1.46 0.25]);
mean_vec_b = vec([4.75 1.46]);
mean_vec_c = vec([5.0 1.86]);
cov_mat_a = [0.029 0.0055; 0.0055 0.0099];
cov_mat_b = [1.097 0.338; 0.3389 0.1189];
cov_mat_c = [0.35 0.22; 0.21 0.17];
d_a = MvNormal(mean_vec_a,cov_mat_a);
d_b = MvNormal(mean_vec_b,cov_mat_b);
d_c = MvNormal(mean_vec_c,cov_mat_c);

a = collect(0:0.05:8);
b = collect(0:0.05:2.5);
pdf_mv = zeros(length(a),length(b));
for i=1:length(a)
    for j=1:length(b)
        pdf_mv[i,j] = maximum([pdf(d_a,[a[i],b[j]]),pdf(d_b,[a[i],b[j]]),pdf(d_c,[a[i],b[j]])]);
    end
end
myplot = plot(layer(x=data_mat_a[:,1],y=data_mat_a[:,2],
Geom.point,Theme(default_color=colorant"red")),layer(x=data_mat_b[:,1],y=data_mat_b[:,2],
Geom.point,Theme(default_color=colorant"blue")),layer(x=data_mat_c[:,1],y=data_mat_c[:,2],
Geom.point,Theme(default_color=colorant"green")),layer(z=pdf_mv,x=a,y=b, Geom.contour(levels=80),
Coord.Cartesian(xmin=0, xmax=8,ymin=0,ymax=2.55),
major_label_font_size=18pt,
minor_label_font_size=14pt,
key_title_font_size = 18pt,
key_label_font_size = 14pt,
major_label_color=colorant"black",
minor_label_color=colorant"black",Guide.xlabel("PetalLength (in cms)"),Guide.ylabel("PetalWidth (in cmas)"))
#draw(PNG("./figs/fitting_multivar_iris.png", 6inch, 3inch,dpi=300), myplot);

# if command:
# if a <b 
#    command
# end