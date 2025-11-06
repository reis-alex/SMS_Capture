classdef SPART_casadi
    properties 
        % variables exchanged between functions
        robot;
        state_vars
        tau
        R0
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
            obj.state_vars.q      = SX.sym('q', obj.robot.n_q,1);
            obj.state_vars.qdot   = SX.sym('qdot', obj.robot.n_q,1);
            obj.tau   = SX.sym('tau',obj.robot.n_q,1);
            obj.state_vars.theta0  = SX.sym('theta0',3,1);
            obj.R0      = SX.sym('R0',3,3);
            obj.state_vars.r0      = SX.sym('r0', 3,1);
            obj.state_vars.r0dot   = SX.sym('r0dot', 3,1);
            obj.state_vars.omega0  = SX.sym('omega_tt', 3,1);

            kin = obj.kinematicsf();
            obj.kinematics.rL = kin.rLf;
            obj.kinematics.RL = kin.RLf;

            obj.dynamics.ddX = obj.ddXf();
            obj.ddX = obj.ddXf();
            diffkin = obj.diffkinematicsf();

            obj.diffkinematics.t0 = diffkin.t0f;
            obj.diffkinematics.tL = diffkin.tLf;
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
            H = obj.H();
            n_q = size1(H.Hm);
            H = [H.H0 H.H0m;
                 H.H0m' H.Hm];

            C = obj.C();
            C = [C.C0 C.C0m;
                C.Cm0 C.Cm];
            f = inv(H)*([zeros(6,1);obj.tau] - C*[obj.state_vars.omega0;obj.state_vars.r0dot;obj.state_vars.qdot]);
            out = casadi.Function('Robot_acceleration',{obj.R0,obj.state_vars.r0,obj.state_vars.omega0,obj.state_vars.q,obj.state_vars.r0dot,obj.state_vars.qdot,obj.tau},{f},{'R0','r0','omega0','q','r0dot','qdot','tau'},{'ddX'});
        end

        function out = diffkinematicsf(obj)
            %Diferential Kinematics
            kin = obj.kinematicsf();
            [Bij,Bi0,P0,pm]     = DiffKinematics_sym(obj.R0,obj.state_vars.r0,kin.rL,kin.e,kin.g,obj.robot);
            [t0,tL]             = Velocities_sym(Bij,Bi0,P0,pm,[obj.state_vars.omega0; obj.state_vars.r0dot],obj.state_vars.qdot,obj.robot);
            out.t0 = t0;
            out.tL = tL;
            out.P0 = P0;
            out.pm = pm;
            out.Bij = Bij;
            out.Bi0 = Bi0;

            out.t0f = casadi.Function('t0',{obj.R0,obj.state_vars.omega0,obj.state_vars.r0dot},{t0},{'R0','omega0','r0dot'},{'t0'});
            out.tLf = casadi.Function('tL',{obj.R0,obj.state_vars.r0,obj.state_vars.omega0,obj.state_vars.q,obj.state_vars.r0dot,obj.state_vars.qdot},{tL},{'R0','r0','omega0','qm','r0dot','qdot'},{'tL'});
        end

        function out = kinematicsf(obj)
            [~,RL,rJ,rL,e,g]     = Kinematics_sym(obj.R0,obj.state_vars.r0,obj.state_vars.q,obj.robot);
            out.RL = RL;
            out.rJ = rJ;
            out.rL = rL;
            out.e = e;
            out.g = g;

%             if exist('idx','var')
%                 out.RLf = casadi.Function('RL',{obj.R0,obj.state_vars.r0,obj.state_vars.q},RL(:,idx));
%                 out.rLf = casadi.Function('rL',{obj.R0,obj.state_vars.r0,obj.state_vars.q},{rL(:,idx)},{'R0','r0','q'},{'rL'});
%             else
                out.RLf = casadi.Function('RL',{obj.R0,obj.state_vars.r0,obj.state_vars.q},RL);
                out.rLf = casadi.Function('rL',{obj.R0,obj.state_vars.r0,obj.state_vars.q},{rL},{'R0','r0','q'},{'rL'});
%             end

        end

        function out = I(obj)
            % Inertias in inertial frames
            kin = obj.kinematicsf();
            [I0,Im]             =   I_I_sym(obj.R0,kin.RL,obj.robot);
            out.I0 = I0;
            out.Im = Im;
            out.Imf = casadi.Function('Im',{obj.R0,obj.state_vars.q},{Im{end}},{'R0','qm'},{'Im'});
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
            out.Hf = casadi.Function('H',{obj.R0,obj.state_vars.r0,obj.state_vars.q},{H},{'R0','base_pos','joint_pos'},{'H'});
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
            out.Cf = casadi.Function('C',{obj.R0,obj.state_vars.r0,obj.state_vars.q,obj.state_vars.omega0,obj.state_vars.r0dot,obj.state_vars.qdot},{C},{'R0','base_pos','joint_pos','omega0','r0dot','joint_vel'},{'C'});
        end

        function out = Jacob(obj,idx)
            diffkin = obj.diffkinematicsf();
            kin = obj.kinematicsf();
            rp = kin.rL(:,idx);
            [J0,Jm] = Jacob_sym(rp,obj.state_vars.r0,kin.rL,diffkin.P0,diffkin.pm,idx,obj.robot);
            out.J0 = J0;
            out.Jm = Jm;
            out.J0f =  casadi.Function('J0',{obj.R0,obj.state_vars.r0,obj.state_vars.q},{J0},{'R0','r0','qm'},{'J0'});
            out.Jmf = casadi.Function('Jm',{obj.R0,obj.state_vars.r0,obj.state_vars.q},{Jm},{'R0','r0','qm'},{'Jm'});
        end

        function out = Jacobdot(obj,idx)
            diffkin = obj.diffkinematicsf();
            kin = obj.kinematicsf();
            rp = kin.rL(:,idx);
            tp = casadi.SX.sym('tp',6,1); %what is this
            [J0dot, Jmdot] = Jacobdot_sym(rp,tp,obj.state_vars.r0,diffkin.t0,kin.rL,diffkin.tL,diffkin.P0,diffkin.pm,idx,obj.robot);
            out.J0dot = J0dot;
            out.Jmdot = Jmdot;
            out.J0dotf = casadi.Function('J0dot',{obj.R0,obj.state_vars.omega0,obj.state_vars.r0,obj.state_vars.q,obj.state_vars.r0dot,tp},{J0dot},{'R0','omega0','r0','qm','r0dot','tp'},{'J0dot'});
            out.Jmdotf = casadi.Function('Jmdot',{obj.R0,obj.state_vars.omega0,obj.state_vars.r0,obj.state_vars.q,obj.state_vars.r0dot,obj.state_vars.qdot,tp},{Jmdot},{'R0','omega0','r0','qm','r0dot','qdot','tp'},{'Jmdot'});
        end

        function out = I_I_dot(obj)
            diffkin = obj.diffkinematicsf();
            kin = obj.kinematicsf();
            [I0_dot,I_I_d] = I_I_dot_sym(obj.R0, diffkin.t0, kin.RL,diffkin.tL,obj.robot);
            out.I0 = I0_dot;
            out.I_I = I_I_d;
        end

        function out = NOC(obj)
            diffkin = obj.diffkinematicsf();
            kin = obj.kinematicsf();
            out = NOC_sym(obj.state_vars.r0,kin.rL,diffkin.P0,diffkin.pm,obj.robot);
        end

        function out = NOC_dot(obj)
            diffkin = obj.diffkinematicsf();
            kin = obj.kinematicsf();
            out = NOCdot_sym(obj.state_vars.r0,diffkin.t0,kin.rL,diffkin.tL,diffkin.P0,diffkin.pm,obj.robot);
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
            [H0_dot, H0m_dot, Hm_dot] = GIM_NOCdot_parsed_sym(N, Ndot, I0, Im, I0_dot, I_I_d, obj.robot);
            out.H0 = H0_dot;
            out.Hom = H0m_dot;
            Hdotf = [H0_dot H0m_dot; H0m_dot' Hm_dot];
            out.Hdot = casadi.Function('Hdot',{obj.R0,obj.state_vars.omega0,obj.state_vars.r0,obj.state_vars.q,obj.state_vars.r0dot,obj.state_vars.qdot},{Hdotf},{'R0','omega0','r0','qm','r0dot','qdot'},{'Hdot'});
            out.H0_dotf = casadi.Function('H0_dot',{obj.R0,obj.state_vars.omega0,obj.state_vars.r0,obj.state_vars.q,obj.state_vars.r0dot,obj.state_vars.qdot},{H0_dot},{'R0','omega0','r0','qm','r0dot','qdot'},{'H0dot'});
            out.H0m_dotf = casadi.Function('H0m_dot',{obj.R0,obj.state_vars.omega0,obj.state_vars.r0,obj.state_vars.q,obj.state_vars.r0dot,obj.state_vars.qdot},{H0m_dot},{'R0','omega0','r0','qm','r0dot','qdot'},{'H0mdot'});
            out.Hm_dotf = casadi.Function('Hm_dot',{obj.R0,obj.state_vars.omega0,obj.state_vars.r0,obj.state_vars.q,obj.state_vars.r0dot,obj.state_vars.qdot},{Hm_dot},{'R0','omega0','r0','qm','r0dot','qdot'},{'Hmdot'});
           
        end


    end


end