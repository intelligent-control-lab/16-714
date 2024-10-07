%% This function returns the A and B matrics for Linear Dynamic Systems
% name: identifies the dynamic model
% dim:  dimention of the state
% mode: CT or DT
% dt: if in DT
% dmode: if in DT

function [A,B] = getAB(name, varargin)
dim = varargin{1};
mode = varargin{2};
if mode == 'DT'
    dt = varargin{3};
    if length(varargin) > 3
        dmode = varargin{4};
    else
        dmode = 'ZOH';
    end
end
switch name
    case 'single_integrator'
        if mode == 'CT'
            A = zeros(dim,dim); B = eye(dim,dim);
        else
            A = eye(dim,dim); B = dt*eye(dim,dim);
        end
    case 'double_integrator'
        n = dim/2;
        if mode == 'CT'
            A = [zeros(n,n) eye(n,n);zeros(n,n) zeros(n,n)]; 
            B = [zeros(n,n);eye(n,n)];
        else
            switch dmode
                case 'ZOH'
                    A = [eye(n,n) dt*eye(n,n);zeros(n,n) eye(n,n)];
                    B = [dt^2/2*eye(n,n);dt*eye(n,n)];
                otherwise
                    A = [eye(n,n) dt*eye(n,n);zeros(n,n) eye(n,n)];
                    B = [0;dt*eye(n,n)];
            end
        end
    otherwise
        A = zeros(dim,dim); B = zeros(dim,dim);
end
end