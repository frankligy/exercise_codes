using Distributions;
using DataFrames;

# define E-Step function
function E_Step(x,mu_SE,mu_VE,mu_VI,sigma,P,mu)
    numerator = P*pdf.(Normal(mu,sigma),x);
    denom = 0.34 .*pdf.(Normal(mu_SE,sigma),x) .+ 0.33 .* pdf.(Normal(mu_VE,sigma),x) .+ 0.33 .* pdf.(Normal(mu_VI,sigma),x);
    post_x = numerator ./denom;
    return post_x
end

# define M-Step function
function M_Step(x,post_x)
    mu = (post_x'*x)./sum(post_x);
    return mu
end

# define Expectation-Maximization process
function EM(x,mu_SE,mu_VE,mu_VI,sigma)
    maxIter = 1000;
    for i=1:maxIter
        print(i,"\n");
        post_x_SE = E_Step(x,mu_SE,mu_VE,mu_VI,sigma,0.34,mu_SE);
        post_x_VE = E_Step(x,mu_SE,mu_VE,mu_VI,sigma,0.33,mu_VE);
        post_x_VI = E_Step(x,mu_SE,mu_VE,mu_VI,sigma,0.33,mu_VI);
        mu_SE_new = M_Step(x,post_x_SE);
        mu_VE_new = M_Step(x,post_x_VE);
        mu_VI_new = M_Step(x,post_x_VI);
        if(abs(mu_SE-mu_SE_new)<0.001 && abs(mu_VE-mu_SE_new)<0.001 && abs(mu_VI-mu_VI_new)< 0.001)
            break;
            end;
        mu_SE = mu_SE_new;
        mu_VE = mu_VE_new;
        mu_VI = mu_VI_new;
    end
    return mu_SE,mu_VE,mu_VI;
end

# main function
using RDatasets;
iris = dataset("datasets","iris");
x = iris[:PetalLength];
mu_SE = 2.0
mu_VE = 3.0
mu_VI = 5.0
sigma = 0.54;
EM(x,mu_SE,mu_VE,mu_VI,sigma)

# plot hte final result
using Gadfly;
using Cairo;
using Fontconfig;
myplot = plot(layer(iris,x=:PetalLength,color=:Species,Geom.histogram,
                    Theme(default_color=colorant"purple")),
        layer(x=0:0.02:8,y=pdf.(Normal(1.47,0.54),0:0.02:8)*10,Geom.line,
                    Theme(default_color=colorant"blue")),
        layer(x=0:0.02:8,y=pdf.(Normal(4.36,0.54),0:0.02:8)*5,Geom.line,
                    Theme(default_color=colorant"red")),
        layer(x=0:0.02:8,y=pdf.(Normal(5.54,0.54),0:0.02:8)*5,Geom.line,
                    Theme(default_color=colorant"green")),
        Guide.xlabel("PetalLength (in cms)"),Guide.ylabel(""), major_label_font_size=18pt,
        minor_label_font_size=14pt,
        key_title_font_size = 18pt,
        key_label_font_size = 14pt,
        major_label_color=colorant"black",
        minor_label_color=colorant"black",Coord.Cartesian(xmin=0, xmax=8));