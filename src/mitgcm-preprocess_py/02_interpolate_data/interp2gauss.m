function [Tp] = interp2gauss(yg,xg,Tg,yp,xp,mask,scale,radius)
%% the yg and xg should be rows (1:n), (1:m) and Tg (m*n)type!!
%% Gaussian interpolation to find Tp, the values of the underlying
%% 2-D function Tg at the points in matrices xp and yp. Matrices xg
%% and yg specify the points at which the data Z is given.

%% --> Initialization;

nxg = length(xg);
nyg = length(yg);

nxp = length(xp);
nyp = length(yp);

Tp = zeros(nxp,nyp);


%% --> Loop;

for ip = 1:nxp
    for jp = 1:nyp
        
        if (mask(ip,jp)~=0)
            
            wp = 0;
            sp = 0;
            
            igp = find( (xg<=xp(ip)+radius) & (xg>=xp(ip)-radius) );
            jgp = find( (yg<=yp(jp)+radius) & (yg>=yp(jp)-radius) );
            
            for ig = igp
                for jg = jgp
                    
                    if ( (Tg(ig,jg)~=0) )
                        
                        wd = exp(-((xp(ip)-xg(ig))^2+(yp(jp)-yg(jg))^2)/scale);
                        wp = wp + wd;
                        sp = sp + wd*Tg(ig,jg);
                    end
                    
                end
            end
            
            if wp == 0
                Tp(ip,jp) = 0;
            else
                Tp(ip,jp) = sp/wp;
            end
            
        end
        
    end
end