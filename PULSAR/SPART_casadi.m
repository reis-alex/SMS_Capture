classdef SPART_casadi
    properties 
        % variables exchanged between functions
        robot;
        state_vars
        tau
        R0
        idx
        % functions to be generated once, then used in simulation
        ddX
        diffkinematics
        kinematics
        dynamics
        Jacobians
    end
    
    methods
        function obj = SPART_casadi(path)
            [robot,~] = urdf2robot_flex_visu(path);
            obj.robot = robot;
            import casadi.*

            % sort dimensions
                obj.idx.r0 = [4:6];
            if robot.n_q > 0
                obj.idx.q = 7:6+obj.robot.n_q;
                obj.idx.omega0 = obj.idx.q(end)+1:obj.idx.q(end)+3;
                obj.idx.r0dot = obj.idx.omega0(end)+1:obj.idx.omega0(end)+3;
                obj.idx.qdot = obj.idx.r0dot(end)+1:obj.idx.r0dot(end)+robot.n_q;
                obj.idx.velocities = obj.idx.omega0(1):obj.idx.qdot(end);
            else
                obj.idx.q = [];
                obj.idx.omega0 = 7:9;
                obj.idx.r0dot = obj.idx.omega0(end)+1:obj.idx.omega0(end)+3;
                obj.idx.qdot = [];
                obj.idx.velocities = obj.idx.omega0(1):obj.idx.r0dot(end);
            end
                obj.idx.positions = 1:obj.idx.omega0(1)-1;

            % create (unique) symbolic variables
            theta0  = SX.sym('theta',3,1);
            r0      = SX.sym('r', 3,1);
            q      = SX.sym('q', obj.robot.n_q,1);
            omega0  = SX.sym('omega', 3,1);
            r0dot   = SX.sym('rdot', 3,1);
            qdot   = SX.sym('qdot', obj.robot.n_q,1);

            obj.state_vars = vertcat(theta0, r0, q, omega0, r0dot, qdot);
            obj.tau   = SX.sym('tau',obj.robot.n_q,1);
            obj.R0      = SX.sym('R0',3,3);

            kin = obj.kinematicsf();
            obj.kinematics.rL = kin.rLf;
            obj.kinematics.RL = kin.RLf;

            obj.dynamics.ddX = obj.ddXf();
            obj.ddX = obj.ddXf();
            diffkin = obj.diffkinematicsf();

            obj.diffkinematics.t0 = diffkin.t0f;
            obj.diffkinematics.tL = diffkin.tLf;
            obj.diffkinematics.P0 = diffkin.P0f;
            obj.diffkinematics.pm = diffkin.pmf;
            Hs = obj.H();
            obj.dynamics.H = Hs.Hf;
            Cs = obj.C();
            obj.dynamics.C = Cs.Cf;

%             Jacob = obj.Jacob();
%             obj.f_J.J0 = Jacob.J0f;
%             obj.f_J.Jm = Jacob.Jmf;

%             Jacobdot = obj.Jacobdot();
%             obj.f_Jdot.J0dot = Jacobdot.J0dotf;
%             obj.f_Jdot.Jmdot = Jacob.Jmdotf;

            Hdot = obj.Hdot();
            obj.dynamics.Hdot = Hdot.Hdot;
        end

        function out = ddXf(obj)

            theta0  = obj.state_vars(1:3);
            r0      = obj.state_vars(4:6);
            q       = obj.state_vars(obj.idx.q);
            omega0  = obj.state_vars(obj.idx.omega0);
            r0dot    = obj.state_vars(obj.idx.r0dot);
            qdot    = obj.state_vars(obj.idx.qdot);

            H = obj.H();
            n_q = size1(H.Hm);
            H = [H.H0 H.H0m;
                 H.H0m' H.Hm];

            C = obj.C();
            C = [C.C0 C.C0m;
                C.Cm0 C.Cm];
            f = inv(H)*([zeros(6,1);obj.tau] - C*[omega0;r0dot;qdot]);
            out = casadi.Function('Robot_acceleration',{obj.R0,obj.state_vars,obj.tau},{f},{'R0','StateVector','tau'},{'ddX'});
        end

        function out = diffkinematicsf(obj)
            theta0  = obj.state_vars(1:3);
            r0      = obj.state_vars(4:6);
            omega0  = obj.state_vars(obj.idx.omega0);
            r0dot    = obj.state_vars(obj.idx.r0dot);
            qdot    = obj.state_vars(obj.idx.qdot);
            q       = obj.state_vars(obj.idx.q);

            %Diferential Kinematics
            kin = obj.kinematicsf();
            [Bij,Bi0,P0,pm]     = DiffKinematics_sym(obj.R0,r0,kin.rL,kin.e,kin.g,obj.robot);
            [t0,tL]             = Velocities_sym(Bij,Bi0,P0,pm,[omega0; r0dot],qdot,obj.robot);
            out.t0 = t0;
            out.tL = tL;
            out.P0 = P0;
            out.pm = pm;
            out.Bij = Bij;
            out.Bi0 = Bi0;

            out.t0f = casadi.Function('t0',{obj.R0,obj.state_vars},{t0},{'R0','StateVector'},{'t0'});
            out.tLf = casadi.Function('tL',{obj.R0,obj.state_vars},{tL},{'R0','StateVector'},{'tL'});
            out.P0f  = casadi.Function('P0',{obj.R0},{P0},{'R0'},{'P0'});
            out.pmf  = casadi.Function('pm',{obj.R0,obj.state_vars},{pm},{'R0','StateVector'},{'pm'});
        end

        function out = kinematicsf(obj)
            theta0  = obj.state_vars(1:3);
            r0      = obj.state_vars(4:6);
            q       = obj.state_vars(obj.idx.q);
            omega0  = obj.state_vars(obj.idx.omega0);
            r0dot    = obj.state_vars(obj.idx.r0dot);
            qdot    = obj.state_vars(obj.idx.qdot);

            [~,RL,rJ,rL,e,g]     = Kinematics_sym(obj.R0,r0,q,obj.robot);
            out.RL = RL;
            out.rJ = rJ;
            out.rL = rL;
            out.e = e;
            out.g = g;

%             if exist('idx','var')
%                 out.RLf = casadi.Function('RL',{obj.R0,obj.state_vars.r0,obj.state_vars.q},RL(:,idx));
%                 out.rLf = casadi.Function('rL',{obj.R0,obj.state_vars.r0,obj.state_vars.q},{rL(:,idx)},{'R0','r0','q'},{'rL'});
%             else
                out.RLf = casadi.Function('RL',{obj.R0,obj.state_vars},RL);
                out.rLf = casadi.Function('rL',{obj.R0,obj.state_vars(1:obj.robot.n_q+6)},{rL},{'R0','StateVector (positions)'},{'rL'});
%             end

        end

        function out = I(obj)
            theta0  = obj.state_vars(1:3);
            r0      = obj.state_vars(4:6);
            omega0  = obj.state_vars(obj.idx.omega0);
            r0dot    = obj.state_vars(obj.idx.r0dot);
            qdot    = obj.state_vars(obj.idx.qdot);
            q       = obj.state_vars(obj.idx.q);

            % Inertias in inertial frames
            kin = obj.kinematicsf();
            [I0,Im]             =   I_I_sym(obj.R0,kin.RL,obj.robot);
            out.I0 = I0;
            out.Im = Im;
            out.Imf = casadi.Function('Im',{obj.R0,q},{Im{end}},{'R0','StateVector'},{'Im'});
        end

        function out = M(obj)
            I = obj.I();
            diffkin = obj.diffkinematicsf();
            %Mass Composite Body matrix
            [M0_tilde,Mm_tilde] =   MCB_sym(I.I0,I.Im,diffkin.Bij,diffkin.Bi0,obj.robot);
            out.M0 = M0_tilde;
            out.Mm = Mm_tilde;
        end

        function out = H(obj)
            M = obj.M();
            diffkin = obj.diffkinematicsf();
            %Generalized Inertia matrix
            [H0,H0m,Hm]         =   GIM_sym(M.M0,M.Mm,diffkin.Bij,diffkin.Bi0,diffkin.P0,diffkin.pm,obj.robot);
            out.H0 = H0;
            out.H0m = H0m;
            out.Hm = Hm;
            H = [H0 H0m;
                 H0m' Hm];
            out.Hf = casadi.Function('H',{obj.R0,obj.state_vars},{H},{'R0','StateVector'},{'H'});
        end
        
        function out = C(obj)
            diffkin = obj.diffkinematicsf();
            I = obj.I();
            M = obj.M();
            %Generalized Convective Inertia matrix
            [C0, C0m, Cm0, Cm]  =   CIM_sym(diffkin.t0,diffkin.tL,I.I0,I.Im,M.M0,M.Mm,diffkin.Bij,diffkin.Bi0,diffkin.P0,diffkin.pm,obj.robot);
            out.C0 = C0;
            out.C0m = C0m;
            out.Cm0 = Cm0;
            out.Cm = Cm;
            C = [C0 C0m; Cm0 Cm];
            out.Cf = casadi.Function('C',{obj.R0,obj.state_vars},{C},{'R0','StateVector'},{'C'});
        end

        function out = Jacob(obj,idx)
            theta0  = obj.state_vars(1:3);
            r0      = obj.state_vars(4:6);
            omega0  = obj.state_vars(obj.idx.omega0);
            r0dot    = obj.state_vars(obj.idx.r0dot);
            qdot    = obj.state_vars(obj.idx.qdot);
            q       = obj.state_vars(obj.idx.q);

            P0 = obj.diffkinematics.P0(obj.R0);
            pm = obj.diffkinematics.pm(obj.R0,obj.state_vars);
            rL = obj.kinematics.rL(obj.R0,obj.state_vars(1:length(obj.state_vars)/2));
            [J0,Jm] = Jacob_sym(rL(:,idx),r0,rL,P0,pm,idx,obj.robot);
            out.J0 = J0;
            out.Jm = Jm;
            out.J0f = casadi.Function('J0',{obj.R0,obj.state_vars},{J0},{'R0','StateVector'},{'J0'});
            out.Jmf = casadi.Function('Jm',{obj.R0,obj.state_vars},{Jm},{'R0','StateVector'},{'Jm'});
        end

        function out = Jacobdot(obj,idx)
            theta0  = obj.state_vars(1:3);
            r0      = obj.state_vars(4:6);
            omega0  = obj.state_vars(obj.idx.omega0);
            r0dot    = obj.state_vars(obj.idx.r0dot);
            qdot    = obj.state_vars(obj.idx.qdot);
            q       = obj.state_vars(obj.idx.q);

            P0 = obj.diffkinematics.P0(obj.R0);
            pm = obj.diffkinematics.pm(obj.R0,obj.state_vars);
            rL = obj.kinematics.rL(obj.R0,obj.state_vars);
            tL = obj.diffkinematics.tL(obj.R0,obj.state_vars);
            t0 = obj.diffkinematics.t0(obj.R0,obj.state_vars);

            [J0dot, Jmdot] = Jacobdot_sym(rL(:,idx),tL(:,idx),r0,t0,rL,tL,P0,pm,idx,obj.robot);
            out.J0dot = J0dot;
            out.Jmdot = Jmdot;
            out.J0dotf = casadi.Function('J0dot',{obj.R0,obj.state_vars},{J0dot},{'R0','StateVector'},{'J0dot'});
            out.Jmdotf = casadi.Function('Jmdot',{obj.R0,obj.state_vars},{Jmdot},{'R0','StateVector'},{'Jmdot'});
        end

        function out = I_I_dot(obj)
            diffkin = obj.diffkinematicsf();
            kin = obj.kinematicsf();
            [I0_dot,I_I_d] = I_I_dot_sym(obj.R0, diffkin.t0, kin.RL,diffkin.tL,obj.robot);
            out.I0 = I0_dot;
            out.I_I = I_I_d;
        end

        function out = NOC(obj)
            r0      = obj.state_vars(4:6);
            diffkin = obj.diffkinematicsf();
            kin = obj.kinematicsf();
            out = NOC_sym(r0,kin.rL,diffkin.P0,diffkin.pm,obj.robot);
        end

        function out = NOC_dot(obj)
            r0      = obj.state_vars(4:6);
            diffkin = obj.diffkinematicsf();
            kin = obj.kinematicsf();
            out = NOCdot_sym(r0,diffkin.t0,kin.rL,diffkin.tL,diffkin.P0,diffkin.pm,obj.robot);
        end

        function out = Hdot(obj)
            N = obj.NOC();
            Ndot = obj.NOC_dot();
            I             =   obj.I();
            I0 = I.I0;
            Im = I.Im;
            Idot = obj.I_I_dot();
            I0_dot = Idot.I0;
            I_I_d = Idot.I_I;
            [H0_dot, H0m_dot, Hm_dot] = GIM_NOCdot_parsed_sym(N, Ndot, I0, Im, I0_dot, I_I_d, obj.robot); %tester avec GIM_sym?
            out.H0 = H0_dot;
            out.Hom = H0m_dot;
            Hdotf = [H0_dot H0m_dot; H0m_dot' Hm_dot];
            out.Hdot = casadi.Function('Hdot',{obj.R0,obj.state_vars},{Hdotf},{'R0','StateVector'},{'Hdot'});
            out.H0_dotf = casadi.Function('H0_dot',{obj.R0,obj.state_vars},{H0_dot},{'R0','StateVector'},{'H0dot'});
            out.H0m_dotf = casadi.Function('H0m_dot',{obj.R0,obj.state_vars},{H0m_dot},{'R0','StateVector'},{'H0mdot'});
            out.Hm_dotf = casadi.Function('Hm_dot',{obj.R0,obj.state_vars},{Hm_dot},{'R0','StateVector'},{'Hmdot'});
           
        end


    end


end