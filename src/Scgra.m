classdef(Abstract) Scgra < handle
    % SCGRA Sequential Conjugate Gradient-Restoration Algorithm
    % SCGRAによる最適制御の数値解法アルゴリズム
    % それぞれの問題に合わせ、このクラスを継承して抽象プロパティ、抽象メソッドを実装すること
    % 初期推定解を用意しておいて、コンストラクタ、Init(), SetInitialGuess(), Calc()の順でメソッドを呼ぶとSCGRAが実行される
    % 具体的な使い方はCar.mやRunCar.mを参照のこと
    
    % References
    % [1] Wu, A. K., and Angelo Miele. "Sequential conjugate gradient-restoration algorithm for optimal control problems with non-differential constraints and general boundary conditions, part I.", Optimal Control Applications and Methods 1.1, 1980, pp.69-88
    % [2] 原田正範、"最適制御問題による地面付近における滑空機の距離最大飛行の解析に関する研究"、東海大学学位論文、1996
    % [3] 麥谷高志、"状態量不等式拘束を伴う最適問題の数値解法に関する研究"、東京大学博士論文、1996
    
    properties
        x;      % n * (N+1)
        u;      % m * (N+1)
        pi;     % p * 1     時間不変を仮定
        lambda; % n * (N+1)
        mu;     % q * 5
        rho;    % k * (N+1)
        
        % Temporary variables
        x_old; % n * (N+1)
        u_old; % m * (N+1)
        pi_old;% p * 1
        
        % 1ステップ前の値
        A_hat; % n * (N+1)  dx/alpha
        B_hat; % m * (N+1)  du/alpha
        C_hat; % p * 1     dpi/alpha
        
        % Gradient Phaseでgammaが計算可能か判断する
        firstRunFlag = 1;
        IterCnt = 0;
    end
    
    properties(Abstract)
        n; % 状態ベクトルの個数
        m; % 入力ベクトルの個数
        p; % パラメータの個数
        k; % 入力拘束の数
        q; % 終端拘束の数
        e1; % epsilon_1 Pの許容値
        e2; % epsilon_2 Qの許容値
        maxIter; %calcの最大繰り返し数
        eps_gs; %黄金分割方法の許容値
        
        Past; % P_* ：Gradient PhaseにおけるPの許容値
        N;%分割数
        
        conjugateFlag;% 0：勾配法、1:共役勾配法
        saveEnableFlag;
    end
    
    methods(Static)
        function ret = Squared(x)
           tmp = x; %2回xが評価されると重くなりそうなので変数においておく
           ret = tmp.' * tmp;
        end
        
        function yint = rk4(F,dt,t,y) % Runge-Kutta 4th
            k_1 = F(t,y);
            k_2 = F(t+0.5*dt,y+0.5*dt*k_1);
            k_3 = F((t+0.5*dt),(y+0.5*dt*k_2));
            k_4 = F((t+dt),(y+k_3*dt));
            
            yint = y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;  % main equation
        end
        
        function fdot = FiniteDiff(f,h)
            n_ = size(f,2);
            fdot = zeros(1,n_);
            fdot(1) = (-3 * f(1) + 4*f(2) - f(3)) / (2*h);
            fdot(n_) = (3 * f(n_) - 4*f(n_-1) + f(n_-2)) / (2*h) ;
%             fdot(1) = (f(2) - f(1)) / h;
%             fdot(n_) = (f(n_) - f(n_-1)) / h ;
            for i = 2: n_-1
                fdot(i) = (f(i+1) - f(i-1))/(2*h);
            end
        end
    end 
    
    methods
        function self = Scgra()
        end
        function Init(self)
            self.firstRunFlag = 1;
            self.IterCnt = 0;
            self.A_hat = zeros(self.n,self.N+1);
            self.B_hat = zeros(self.m,self.N+1);
            self.C_hat = zeros(self.p,1);
        end
        
        function SetNum(self,n_in,m_in,p_in,k_in,q_in)
            if n_in <= 0 || m_in <= 0 || p_in <= 0 || k_in <= 0 || q_in <= 0
                error('Init error');
            end
  
            self.n = n_in;
            self.m = m_in;
            self.p = p_in;
            self.k = k_in;
            self.q = q_in;
        end
        
        function SetDivisionNumber(self,N_in)
            if N_in <= 0
                error('SetDivisionNumber error');
            end
            self.N = N_in;
        end
        
        function SetInitialGuess(self,x_est, u_est, pi_est)
            if size(x_est,1) ~= self.n || size(x_est,2) ~= self.N+1 || size(u_est,1) ~= self.m || size(u_est,2) ~= self.N+1 ||size(pi_est,1) ~= self.p || size(pi_est,2) ~= 1
                error('input size error');
            end
            self.x = x_est;
            self.u = u_est;
            self.pi = pi_est;
        end
        
        function SetConvParam(self, e1_in, e2_in)
            if e1_in <= 0 || e2_in <= 0
                 error('SetConvParam error');
            end
            self.e1 = e1_in;
            self.e2 = e2_in;
        end
        
        function SetMaxIter(self, maxIter_in)
            if maxIter_in <= 0
                error('SetMaxIter error');
            end
            self.maxIter = maxIter_in;
        end
        
        function J = AugEvalFunc(self,x,u,pi,lambda,mu,rho)
            [~,tmp1] = ode45(@(t,dummy) self.AugEvalFuncDiff(t,x,u,pi,lambda,rho), [0 1], 0);
            J = tmp1(end) + mu.' * self.Psi(x(:,end),pi) + self.g(x(:,end),pi);
        end
        
        function fq = Interp(self,t,f)
            % f: n * (N+1)
            n_ = size(f,1);
            fq = zeros(n_,1);
            N_ = self.N;
            tary = linspace(0,1,N_+1);
            for i=1:n_
                fq(i,1) = interp1(tary, f(i,:), t);
            end
        end
        
        function fq = InterpDiff(self,t,f)
            % f: n * (N+1)
            n_ = size(f,1);
            fq = zeros(n_,1);
%             fdot = diff(f,1,2) * self.N; % n * N  引数は要チェック 列間の一次微分をとりたい
%             stride = 1/self.N /2;
            N_ = self.N;
            tary = linspace(0,1,N_+1);
            for i=1:n_
                fdot = Scgra.FiniteDiff(f(i,:),1/N_);
                fq(i,1) = interp1(tary,fdot,t,'linear','extrap');
            end
        end
                
        function ret = AugEvalFuncDiff(self,t,x,u,pi,lambda,rho)
            xq = self.Interp(t,x);
            uq = self.Interp(t,u);
            lambdaq = self.Interp(t,lambda);
            rhoq = self.Interp(t,rho);
            xdotq = self.InterpDiff(t,x);
            
            ret = self.f(t,xq,uq,pi) + lambdaq.' * (xdotq-self.Phi(t,xq,uq,pi)) + rhoq.' * self.S(t,xq,uq,pi);
        end
        
        function P = ConErr(self)
            P = self.ConErr_(self.x,self.u,self.pi);
        end
        
        function P = ConErr_(self,x,u,pi) %Constraint error P
            % Refer [3]
            [~,tmp1] = ode45(@(t,dummy) Scgra.Squared(self.StateEqErrDiff(t,x,u,pi)), [0 1], 0);
            
            [~,tmp2] = ode45(@(t,dummy) Scgra.Squared(self.EqualityConstraintErr(t,x,u,pi)), [0 1], 0);
            
            tmp3 = self.Psi(x(:,end),pi);
            
            P = tmp1(end) + tmp2(end) + dot(tmp3,tmp3);
            if(P < 0) %数値誤差でPが負となることがあったので念の為に正にしておく
               P = 0; 
            end
        end
        
        function ret = EqualityConstraintErr(self,t,x,u,pi)
            xq = self.Interp(t,x);
            uq = self.Interp(t,u);
            ret = self.S(t,xq,uq,pi);
        end
        
        function ret = StateEqErrDiff(self,t,x,u,pi) %状態方程式のエラー
%             xq = self.Interp(t,x);
%             uq = self.Interp(t,u);
%             xdotq = self.InterpDiff(t,x);
%             ret = xdotq - self.Phi(t,xq,uq,pi);
            xq = self.Interp(t,x);
            if t > 0
                [~,tmp] = ode45(@(t2,dum) self.StateEqInterp(t2, x, u, pi), [0 t], x(:,1));
                ret = xq - tmp(end,:).';
            else
                ret = xq - x(:,1);
            end
        end
        
        function ret = StateEqInterp(self,t2,x,u,pi)
            xq = self.Interp(t2,x);
            uq = self.Interp(t2,u);
            ret = self.Phi(t2,xq,uq,pi);
        end
        
        function Q = OptErr(self) %Error in optimality condition
            [~,tmp1] = ode45(@(t,dummy) Scgra.Squared(self.HxErrDiff(t,self.x,self.u,self.pi,self.lambda,self.rho)), [0 1], 0);
            [~,tmp2] = ode45(@(t,dummy) Scgra.Squared(self.HuInterp(t,self.x,self.u,self.pi,self.lambda,self.rho)'), [0 1], 0);
            if(self.p > 0)
                [~,tmp] = ode45(@(t,dummy) self.HpiInterp(t,self.x,self.u,self.pi,self.lambda,self.rho).', [0 1], zeros(self.p,1));
                tmp3 = tmp(end,:).' + self.Psipi(self.x(:,end),self.pi).'*self.mu + self.gpi(self.x(:,end),self.pi);
            else
                tmp3 = 0;
            end
            tmp4 = self.lambda(:,end) + self.Psix(self.x(:,end),self.pi).'*self.mu + self.gx(self.x(:,end),self.pi)';
            
            Q = tmp1(end) + tmp2(end) + dot(tmp3,tmp3) + dot(tmp4,tmp4);
            if(Q < 0) %数値誤差でPが負となることがあったのでQも念の為に正にしておく
               Q = 0; 
            end  
        end
        
        function ret = HpiInterp(self,t2,x,u,pi,lambda,rho)
            xq = self.Interp(t2,x);
            uq = self.Interp(t2,u);
            lambdaq = self.Interp(t2,lambda);
            rhoq = self.Interp(t2,rho);
            
            ret = self.Hpi(t2,xq,uq,pi,lambdaq,rhoq);
        end
        
        function ret = HuInterp(self,t2,x,u,pi,lambda,rho)
            xq = self.Interp(t2,x);
            uq = self.Interp(t2,u);
            lambdaq = self.Interp(t2,lambda);
            rhoq = self.Interp(t2,rho);
            
            ret = self.Hu(t2,xq,uq,pi,lambdaq,rhoq);
        end
        
        function ret = HxErrDiff(self,t,x,u,pi, lambda, rho)
            lambdaq = self.Interp(t,lambda);
            if t > 0
                [~,tmp] = ode45(@(t2,dum) self.HxInterp(t2,x,u,pi,lambda,rho).', [0 t],lambda(:,1));
                ret = lambdaq - tmp(end,:).';
            else
                ret = lambdaq - lambda(:,1);
            end
        end
        
        function ret = HxInterp(self,t2,x,u,pi,lambda,rho)
            xq = self.Interp(t2,x);
            uq = self.Interp(t2,u);
             lambdaq = self.Interp(t2,lambda);
            rhoq = self.Interp(t2,rho);
            
            ret = self.Hx(t2,xq,uq,pi,lambdaq,rhoq);
        end
        
        function Push(self, x_in, u_in, pi_in)
            self.x_old = x_in;
            self.u_old = u_in;
            self.pi_old = pi_in;
        end
        
        function [x, u, pi] = Pop(self)
            x = self.x_old;
            u = self.u_old;
            pi = self.pi_old;
        end
        
        function Q = OptErrGrad(self,B,C)
            [~,tmp] = ode45(@(t,dummy) self.OptErrGradDiff(t,B), [0 1], 0);
            Q = tmp(end) + C.' * C;
        end
        
        function ret = OptErrGradDiff(self,t,B)
            Bq = self.Interp(t,B);
            ret = Bq.' * Bq;
        end
        
        function Z = calcZ(self,B,C,B_hat,C_hat)
            [~,tmp] = ode45(@(t,dummy) self.ZDiff(t,B,B_hat), [0 1], 0);
            Z = tmp(end) + C_hat.' * C;
        end
        
        function ret = ZDiff(self,t,B,B_hat)
            Bq = self.Interp(t,B);
            B_hatq = self.Interp(t,B_hat);
            
            ret = B_hatq.' * Bq;
        end
        
        function [A, B, C, lambda, mu, rho, alpha] = ConjGrad(self)
            [A,B,C,lambda,mu,rho] = self.SolveLTPBVP(1,0);
            if(self.firstRunFlag == 0)
                Q = self.OptErrGrad(B,C);
                gamma = Q / self.OptErrGrad(self.B_hat,self.C_hat);
                Z = self.calcZ(B,C,self.B_hat,self.C_hat);
                Jalpha0 = -(Q + gamma * Z);
                if(Jalpha0 >= 0)
                    gamma = 0;
                end
            else %First run
                gamma =  0;
                if(self.conjugateFlag == 1)
                    self.firstRunFlag = 0;
                end
            end
            
            A = A + gamma * self.A_hat;
            B = B + gamma * self.B_hat;
            C = C + gamma * self.C_hat;
            
            self.A_hat = A;
            self.B_hat = B;
            self.C_hat = C;
            
            alpha = self.FindOptAlpha_Con(A,B,C,lambda,mu,rho);
            self.Update(alpha, A,B,C,lambda,mu,rho);
            if self.p > 0
                fprintf('Finish Gradient Phase, tau:%f\n', self.pi(1));
            else
                fprintf('Finish Gradient Phase\n');
            end
        end
        
        function [A, B, C, lambda, mu, rho, alpha] = Restore(self)
%             P = self.ConErr();
%             while P > self.e1
%                 [A,B,C,lambda,mu,rho] = self.SolveLTPBVP(0,1);
%                 self.Update(alpha, A,B,C,lambda,mu,rho);
%                 alpha = alpha/2;%self.FindOptAlpha_Res(A, B, C);
%                 P = self.ConErr();
%             end
            A = zeros(self.n,self.N+1);
            B = zeros(self.m,self.N+1);
            C = zeros(self.p,1); 
            lambda  = zeros(self.n,self.N+1);
            mu  = zeros(self.q,1);
            rho  = zeros(self.k,self.N+1);
            alpha = 1;
                
            P = self.ConErr();
            cnt = 0;
            while P > self.e1 && cnt < 10
                P0 = P;
                [A,B,C,lambda,mu,rho] = self.SolveLTPBVP(0,1);
                alpha = 1;
                P = self.ConErr_(self.x+alpha*A,self.u+alpha*B,self.pi+alpha*C);
                cnt2 = 0;
                while P > P0 && cnt2 < 100
                    alpha = alpha/2;
                    P = self.ConErr_(self.x+alpha*A,self.u+alpha*B,self.pi+alpha*C);
                    cnt2 = cnt2+1;
                    fprintf('Rest Cnt2: %d, P:%f, P0:%f \n', cnt2, P,P0);
                end
                if(cnt2 == 0)
                    fprintf('Rest Cnt2: %d, P:%f, P0:%f \n', cnt2, P,P0);
                end
                if(cnt2==15)
                    error('Restore Error: P is not smaller than P0\nP:%f\n', P);
%                     fprintf('Restore Error: P is not smaller than P0\nP:%f\n', P);
                end

                self.Update(alpha, A,B,C,lambda,mu,rho);
            end
            cnt = cnt + 1;
            if(cnt==10)
               error('Restore Error: P cannnot be made small.\nP:%f\n', P);
            end
           if(self.p > 0 ) 
               fprintf('Finish Restoration Phase, tau:%f\n', self.pi(1));
           else 
               fprintf('Finish Restoration Phase\n');
           end
        end
        
        function [Ai, Bi, Ci, lambdai, rhoi, tmp] = paralellFunc(self, wi, kg, kr)
            n_ = self.n;
            m_ = self.m;
            k_ = self.k;
            p_ = self.p;
            q_ = self.q;
            N_ = self.N;
            
            Ai = zeros(n_,N_+1);
            Bi = zeros(m_,N_+1);
            lambdai = zeros(n_,N_+1);
            rhoi = zeros(k_,N_+1);
            tmp = zeros(n_+p_+1+q_,1);            
            
            Ci = wi(n_+1:end,1);
            lambdai(:,1) = wi(1:n_,1);
            
            ALambdai = zeros(2*n_,N_+1); %Aとλは一括で計算する
            ALambdai(:,1) = [Ai(:,1);lambdai(:,1)];
            
            for ts = 0:N_
                t = ts/N_;
                dt = 1/N_;
                
                xq = self.x(:,ts+1);
                uq = self.u(:,ts+1);
                     
                Ai(:,ts+1) = ALambdai(1:n_,ts+1);
                lambdai(:,ts+1) = ALambdai(n_+1:end,ts+1);
                
                [Bi(:,ts+1),rhoi(:,ts+1)] = self.Brho(ts,xq,uq,self.pi,Ai(:,ts+1),Ci,lambdai(:,ts+1),kg,kr);
                
                if(ts == N_)
                    break;
                end
                
                ALambdai(:,ts+2) = Scgra.rk4(@(t2,ALambda) self.ALambdaDiff(t2, self.x, self.u, self.pi, ALambda(1:n_),Ci(:),ALambda(n_+1:2*n_),kg,kr), dt, t, ALambdai(:,ts+1));
            end
            
            tmp(1:q_) = self.Psix(self.x(:,end),self.pi) * Ai(:,end) + self.Psipi(self.x(:,end),self.pi) * Ci;
            if(p_ > 0)
                [~,integral] = ode45(@(t2,dummy) self.tmpDiff(t2,self.x,self.u,self.pi,lambdai(:,:),rhoi(:,:)) ,[0 1],zeros(p_,1));
                tmp(q_+1:q_+p_) = integral(end,:)';
                tmp(q_+1:q_+p_) = tmp(q_+1:q_+p_) + Ci;
            end
            tmp(q_+p_+1:end-1) = lambdai(:,end);
        end
        
        function [A,B,C,lambda,mu,rho] = SolveLTPBVP(self,kg,kr)
            % Linear Two-Point Boundary Value Problems
            n_ = self.n;
            m_ = self.m;
            k_ = self.k;
            p_ = self.p;
            q_ = self.q;
            N_ = self.N;
            
            tmp = zeros(n_+p_+1+q_);
            Ai = zeros(n_,N_+1,n_+p_+1);
            Bi = zeros(m_,N_+1,n_+p_+1);
            Ci = zeros(p_,n_+p_+1);
            lambdai = zeros(n_,N_+1,n_+p_+1);
            rhoi = zeros(k_,N_+1,n_+p_+1);
            
            parfor i=1:n_+p_+1       
                wi=zeros(n_+p_,1);
                if i < n_+p_+1
                    wi(i,1) = 1;
                end
                [Ai(:,:,i), Bi(:,:,i), Ci(:,i), lambdai(:,:,i), rhoi(:,:,i), tmp(:,i)] = self.paralellFunc(wi, kg, kr);
            end
            
            tmp(q_+1:q_+p_,n_+p_+1+1:end) = self.Psipi(self.x(:,end),self.pi).';
            tmp(q_+p_+1:q_+p_+n_,n_+p_+1+1:end) = self.Psix(self.x(:,end),self.pi).';
            tmp(end,1:n_+p_+1) = ones(1,n_+p_+1);
            tmp2 = zeros(n_+p_+1+q_,1);
            tmp2(1:q_,1) = - kr * self.Psi(self.x(:,end),self.pi);
            if kg == 1 && p_ > 0
                [~,integral] = ode45(@(t,dummy) self.fpiInterp(t,self.x,self.u,self.pi,kg) ,[0 1],zeros(p_,1));
                tmp2(q_+1:q_+p_,1) = integral(end,:)' - kg * self.gpi(self.x(:,end),self.pi).';
            elseif p_ > 0
                tmp2(q_+1:q_+p_,1) = zeros(self.p,1);
            end
            tmp2(q_+p_+1:end-1) = -self.gx(self.x(:,end),self.pi).';
            tmp2(end,1) = 1;
            kmu = tmp \ tmp2; %inv(tmp) * tmp2;

            kcoef = kmu(1:n_+p_+1,1);
            mu = kmu(n_+p_+1+1:end,1);
            

            C = Ci * kcoef;
            A = zeros(n_,N_+1);
            B = zeros(m_,N_+1);
            lambda = zeros(n_,N_+1);
            rho = zeros(k_,N_+1);
            
            parfor ts = 1:N_+1
                A(:,ts) = reshape(Ai(:,ts,:),[n_ n_+p_+1]) * kcoef;
                B(:,ts) = reshape(Bi(:,ts,:),[m_ n_+p_+1]) * kcoef;
                lambda(:,ts) = reshape(lambdai(:,ts,:),[n_ n_+p_+1]) * kcoef;
                rho(:,ts) = reshape(rhoi(:,ts,:),[k_ n_+p_+1]) * kcoef;
            end            
        end
        
        function ret = tmpDiff(self,t,x,u,pi,lambda,rho) % Hpi Interpolation Transposed
            xq = self.Interp(t,x);
            uq = self.Interp(t,u);
            lambdaq = self.Interp(t,lambda);
            rhoq = self.Interp(t,rho);
            ret = - self.Phipi(t,xq,uq,pi).'*lambdaq + self.Spi(t,xq,uq,pi).'*rhoq;
        end
        
        function ret = fpiInterp(self,t,x,u,pi,kg)
            xq = self.Interp(t,x);
            uq = self.Interp(t,u);
            ret = - kg * self.fpi(t,xq,uq,pi).';
        end
        
        function [B,rho] = Brho(self,t,x,u,pi,A,C,lambda,kg,kr)
            % x,u,piなどは時刻tにおける値とし、interpはしない
            % x: n * 1, u: m * 1, pi: p * 1
            m_ = self.m;
            k_ = self.k;
            
            left = zeros(m_+k_, m_+k_);
            right = zeros(m_+k_,1);
            left(1:m_,1:m_) = eye(m_,m_);
            left(m_+1:end,m_+1:end) = zeros(k_,k_);

            left(1:m_,m_+1:end) = self.Su(t,x,u,pi).';
            left(m_+1:end,1:m_) = left(1:m_,m_+1:end).';
            right(1:m_) = -kg * self.fu(t,x,u,pi).' + self.Phiu(t,x,u,self.pi).'*lambda;
            right(m_+1:end) = -(kr * self.S(t,x,u,pi) + self.Sx(t,x,u,pi) * A + self.Spi(t,x,u,pi) * C);
            
            Brho = left \ right;%inv(left) * right
            B = Brho(1:m_);
            rho = Brho(m_+1:end,1);
        end
        
        function ret = ALambdaDiff(self,t,x,u,pi,A,C,lambda,kg,kr)
            xq = self.Interp(t,x);
            uq = self.Interp(t,u);
            xdotq = self.InterpDiff(t,x);
            %Aとλはそれぞれ時刻tにおけるn*1行列なのでinterpしない
            
            [B,rho_] = self.Brho(t,xq,uq,pi,A,C,lambda,kg,kr);
            
            ADot =self.Phix(t,xq,uq,pi)*A + self.Phiu(t,xq,uq,pi) * B + self.Phipi(t,xq,uq,pi) * C - kr * (xdotq - self.Phi(t,xq,uq,pi));
            lambdaDot = kg * self.fx(t,xq,uq,pi).' - self.Phix(t,xq,uq,pi).'*lambda + self.Sx(t,xq,uq,pi).'*rho_;
            
            ret = [ADot;lambdaDot];
        end
        
        function alpha = FindOptAlpha_Con(self, A, B, C, lambda, mu, rho)
            alpha = self.GoldenSelection(0,1,@(alpha) self.AugEvalFunc(self.x+alpha*A,self.u+alpha*B,self.pi+alpha*C,lambda,mu,rho));
            P = self.ConErr_(self.x+alpha*A,self.u+alpha*B,self.pi+alpha*C);
%             fprintf('Grad P:%f\n',P);
            while(P > self.Past || (self.p > 0 && self.pi(1)+alpha*C <= 0))
                alpha = alpha/2;
                P = self.ConErr_(self.x+alpha*A,self.u+alpha*B,self.pi+alpha*C);
%                 fprintf('Grad P:%f\n',P);
            end
        end
        
        function alpha = FindOptAlpha_Res(self, A, B, C)
            alpha = self.GoldenSelection(0,1,@(alpha) self.ConErr_(self.x+alpha*A,self.u+alpha*B,self.pi+alpha*C));
        end
        
        function Update(self, alpha, A, B, C, lambda, mu, rho)
            self.x = self.x + alpha * A;
            self.u = self.u + alpha * B;
            self.pi = self.pi + alpha * C;
            
            %lambda, mu, rho をメンバ変数へ保存
            self.lambda = lambda;
            self.mu = mu;
            self.rho = rho;
        end
        
        function ret = H(self,t,x,u,pi,lambda,rho) %Hamiltonian
            ret = self.f(t,x,u,pi) - lambda.' * self.Phi(t,x,u,pi) + rho.' * self.S(t,x,u,pi);
        end
        
        function  ret = Hx(self,t,x,u,pi,lambda,rho)
            ret = self.fx(t,x,u,pi) - lambda.' * self.Phix(t,x,u,pi) + rho.' * self.Sx(t,x,u,pi);
        end
        
        function  ret = Hu(self,t,x,u,pi,lambda,rho)
            ret = self.fu(t,x,u,pi) - lambda.' * self.Phiu(t,x,u,pi) + rho.' * self.Su(t,x,u,pi);
        end
        
        function  ret = Hpi(self,t,x,u,pi,lambda,rho)
            ret = self.fpi(t,x,u,pi) - lambda.' * self.Phipi(t,x,u,pi) + rho.' * self.Spi(t,x,u,pi);
        end
        
        function Calc(self)
            
            P = self.e1 +1;
            Q = self.e2 +1;
            
%             P = self.ConErr();
%             fprintf('Initial P:%f\n',P);

            % Restoration Phase
            [~, ~, ~, ~, ~, ~, ~] = self.Restore();
            
            while (P > self.e1 || Q > self.e2) && (self.IterCnt < self.maxIter)
                I1 = self.EvalFunc(self.x,self.u,self.pi);
                
                self.Push(self.x, self.u, self.pi);
                
                % Conjugate Gradient Phase
                [A_con, B_con, C_con, lambda_con, mu_con, rho_con, alpha] = self.ConjGrad();
                
%                 P = self.ConErr();
%                 Q = self.OptErr();
%                 fprintf('tf:%.3f,P:%f, Q:%f\n',self.pi(1),P,Q);
                
                % Restoration Phase
                [~, ~, ~, ~, ~, ~, ~] = self.Restore();
                
                I3 = self.EvalFunc(self.x,self.u,self.pi);
                cnt2 = 0;
                while I3 > I1 && cnt2 < 100 %もし評価関数が大きくなってしまったらGradient Phaseのステップサイズを半分にしてUpdateし、Restorationをやり直す
                    if(cnt2>15)
                       [self.x, self.u, self.pi] = self.Pop();
                       error('Error:Restoration Phase makes EvalFunc large');
                       fprintf('Break\n',I1,I3);
                       break; 
                    end
                    fprintf('Restoration Phase makes EvalFunc large, I1:%f, I3:%f\n',I1,I3);
                    [self.x, self.u, self.pi] = self.Pop();
                    alpha = alpha/2;
                    self.Update(alpha, A_con, B_con, C_con, lambda_con, mu_con, rho_con);
                    
                    % Restoration Phase
                    [~, ~, ~, ~, ~, ~, ~] = self.Restore();
                    I3 = self.EvalFunc(self.x,self.u,self.pi);
                    cnt2 = cnt2 + 1;
                end
                
                self.IterCnt = self.IterCnt+1;
                P = self.ConErr();
                Q = self.OptErr();
                if(self.saveEnableFlag==1)
                    filename = self.GenFileName();
                    save(filename, 'self');
                end
                if self.p > 0
                    fprintf('tf:%.3f,P:%f, Q:%f\n',self.pi(1),P,Q);
                else 
                    fprintf('I:%.3f,P:%f, Q:%f\n',I3,P,Q);
                end
            end
            if(self.IterCnt==self.maxIter)
               error('Exceed maxIter'); 
            end
        end
        
        function I = EvalFunc(self,x,u,pi)
            [~,tmp] = ode45(@(t,dummy) self.EvalFuncDiff(t,x,u,pi), [0 1], 0);
            I = tmp(end) + self.g(x(:,end),pi);
        end
        
        function ret = EvalFuncDiff(self,t,x,u,pi)
            xq = self.Interp(t,x);
            uq = self.Interp(t,u);
            ret = self.f(t,xq,uq,pi);
        end
        
        function [min] = GoldenSelection(self,low,up,J)
            %黄金分割法
            r = 0.618033988749895;%(-1+sqrt(5))/2;
            
            cnt = 0;
            while(up-low > self.eps_gs && cnt<1000)
                a1 = low + (1-r) * (up - low);
                a2 = low +    r  * (up-low);
                
                J1 = J(a1);
                J2 = J(a2);
%                 fprintf('J1:%f, J2:%f\n',J1,J2);
                if(J1<J2)
                    up = a2;
                    a2=a1;
                else
                    low = a1;
                    a1 = a2;
                end
                %fprintf(.'%f,%f\r\n.',low,up);
                cnt = cnt+1;
            end
            if(cnt>1000)
                error('error in golden selection');
            end
            min = (low+up)/2;
%             fprintf('Ret:%f\n',min);
        end
    end
    
    methods (Abstract) % Need to implement (Abstract method)
        % 以下の関数のx,u,piは時刻tにおける値でサイズは n*1, m*1, p*1
        
        % 評価関数の被積分項
        f_ = f(self,t,x,u,pi)
        fx_ = fx(self,t,x,u,pi)
        fu_ = fu(self,t,x,u,pi)
        fpi_ = fpi(self,t,x,u,pi)
        
        % 終端コスト
        g_ = g(self,x,pi)
        gx_ = gx(self,x,pi)
        gpi_ = gpi(self,x,pi)
        
        % State Equation
        phi = Phi(self,t,x,u,pi) 
        phix= Phix(self,t,x,u,pi)
        phiu = Phiu(self,t,x,u,pi)
        phipi = Phipi(self,t,x,u,pi)
        
        % Input State Constraint 
        s = S(self,t,x,u,pi) 
        sx = Sx(self,t,x,u,pi)
        su = Su(self,t,x,u,pi)
        spi = Spi(self,t,x,u,pi)
        
        % Terminal Constraint
        psi = Psi(self,x,pi)
        psix = Psix(self,x,pi)           
        psipi = Psipi(self,x,pi)
        
        % Generate Filename
        filename = GenFileName(self)
    end
end
