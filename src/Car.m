classdef Car < Scgra
    %CAROPTIM 車の急加速・急停車問題
    
    properties
        n = 2;
        m = 2;
        p = 1;
        k = 1;
        q = 2;
        e1 = 1.4e-4;
        e2 = 1e-3;  
        maxIter = 1000;
        eps_gs = 1e-4; %黄金分割方法の許容値
        
        Past = 0.01;
        N = 50;
        
        xf  = 1;
        umax = 1;
        
        conjugateFlag = 1;
        saveEnableFlag = 0;
    end
    
    methods
        % 以下の関数のx,u,piは時刻tにおける値でサイズは n*1, m*1, p*1
        
        % 評価関数の被積分項
        function f_ = f(self,t,x,u,pi)
            f_ = pi(1);
        end
        function fx_ = fx(self,t,x,u,pi)
            fx_ = zeros(1,size(x,1));
        end
        function fu_ = fu(self,t,x,u,pi)
            fu_ = zeros(1,size(u,1));
        end
        function fpi_ = fpi(self,t,x,u,pi)
            fpi_ = zeros(1,size(pi,1));
            fpi_(1,1) = 1;  
        end

        % 終端コスト
        function g_ = g(self,x,pi)
            g_ = 0;
        end
        function gx_ = gx(self,x,pi)
            gx_ = zeros(1,size(x,1));
        end
        function gpi_ = gpi(self,x,pi)
            gpi_ = zeros(1,size(pi,1));
        end
        
        % State Equation
        function phi = Phi(self,t,x,u,pi)
            phi = zeros(size(x));
            phi(1,1) = pi(1) * x(2);
            phi(2,1) = pi(1) * u(1);
        end
        function phix= Phix(self,t,x,u,pi)
            phix = zeros(size(x,1),size(x,1));
            phix(1,2) = pi(1);
        end
        function phiu = Phiu(self,t,x,u,pi)
            phiu = zeros(size(x,1),size(u,1));
            phiu(2,1) = pi(1);
        end
        function phipi = Phipi(self,t,x,u,pi)
            phipi = zeros(size(x,1),size(pi,1));
            phipi(1,1) = x(2);
            phipi(2,1) = u(1);
        end
        
        % Input State Constraint
        function s = S(self,t,x,u,pi)
            s = zeros(self.k,1);
            s(1,1) = u(1)-self.umax*(sin(u(2)));
        end
        function sx = Sx(self,t,x,u,pi)
            sx = zeros(self.k,size(x,1));
        end
        function su = Su(self,t,x,u,pi)
            su = zeros(self.k,size(u,1));
            su(1,1) = 1;
            su(1,2) = -self.umax * cos(u(2));
        end
        function spi = Spi(self,t,x,u,pi)
            spi = zeros(self.k,size(pi,1));
        end
        
        % Terminal Constraint
        function psi = Psi(self,x,pi)
            psi = zeros(self.q,1);
            psi(1,1) = x(1) - self.xf;
            psi(2,1) = x(2) - 0;
        end
        function psix = Psix(self,x,pi)
            psix = zeros(self.q,size(x,1));
            psix(1,1) = 1;
            psix(2,1) = 0;
            psix(1,2) = 0;
            psix(2,2) = 1;
        end
        function psipi = Psipi(self,x,pi)
            psipi = zeros(self.q,size(pi,1));
            psipi(1,1) = 0;
        end
        
        % Generate Filename
        function filename = GenFileName(self)
            filename = sprintf('./result/Iter_%03.0f_tf_%.0f.mat', self.IterCnt,self.pi);
        end
    end
end