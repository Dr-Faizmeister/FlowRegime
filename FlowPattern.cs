using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlowRegime
{
    class FlowPattern : IHydrodynamicModel
    {
        public Boolean Linear;
        private Byte FlowRegime; // режим течения
                                 /*
                            Дисперсный             = 1
                            Расслоенный сглаженный = 2
                            Расслоенный волновой = 3
                            Кольцевой-дисперсный   = 4
                            Пузырьковый            = 5
                            Пробковый              = 6
                            Снарядный              = 7
                            Пенный                 = 8
                            */
        double D,                      //Диаметр потока
        rho_l, rho_g,           //Плотности жидкости и газа
        mu_l, mu_g,             //Вязкости жидкости и газа
        theta,                  //Угол наклона скважины (вертикальная - 90, горизонтальная - 0)
        sigma,                  //Коэффициент межфазного поверхностного натяжения
        h_pr,                   //Величина неровности трубы

        A, S,                   // Площадь(A) и периметр(S) потока

        //---------------------
        eps,                    //Точность расчета для расслоенного и кольцевого режимов течения
        g,                      //Ускорение свободного падения
        sd,                     //Константа, необходимая для расчета в механистической модели

        //---------------------
        V_sl, V_sg, V_m,        //Приведенные скорости жидкости, газа и смеси

        E_l_s,                  //Фазовое содержание жидкости в пробке жидкости
        FE,                     //Доля жидкости, находящаяся в виде капель в газовом ядре (кольцевой-дисперсный режим)
        rho_c, mu_c,            //Плотность и вязкость газового ядра
        delta_t,                //Критическое значение безразмерной толщины плёнки (кольцевой-дисперсный режим)

        //---------------------
        E_l, dP_dx,             //Фазовое содержание жидкости, градиент давления

        E_l_1,                  //Фазовое содержание жидкости (дисперсный режим)
        E_l_2_3, y_2_3,         //Фазовое содержание жидкости и угол подъема жидкости (расслоенный режим)
        E_l_4, delta_4,         //Фазовое содержание жидкости и толщина пленки жидкости (кольцевой-дисперсный режим)
        E_l_5,                  //Фазовое содержание жидкости (пузырьковый режим)
        E_l_6_7,                //Фазовое содержание жидкости (пробковый и снарядный режимы течения)
        E_l_8;          //Фазовое содержание жидкости (пенный режим)
                        //---------------------
        double cos_theta, sin_theta; // косинус и синус угла наклона скважины
        private float visc_l;
        private float visc_g;

        public FlowPattern(double p1, double p2, double p3, double p4, double p5, double p6, double p7)
        //Конструктор (с заданием параметров и свойств для расчетов)
        /*
                                                                          p1 - диаметр потока
                                                                          p2 - плотность жидкости
                                                                          p3 - плотность газа
                                                                          p4 - вязкость жидкости
                                                                          p5 - вязкость газа
                                                                          p6 - коэффициент поверхностного натяжения
                                                                          p7 - абсолютная неровность трубы*/
        {
            Boolean k1, k2;
            D = p1;
            rho_l = p2;
            rho_g = p3;
            mu_l = p4;
            mu_g = p5;
            sigma = p6;
            h_pr = p7;

            A = 0.25 * Math.PI * D * D;
            S = Math.PI * D;

            eps = 1e-8;
            g = 9.81;
            sd = 0.06;

            //V_sl := 0.0; //0.01524;0.8524
            V_sl = 0.5;
            //V_sg := 0.0; //10;
            V_sg = 5.0;
            V_m = 0;

            E_l_s = 0;
            FE = 0;
            rho_c = 0;
            mu_c = 0;
            delta_t = 0;

            E_l = 0;
            dP_dx = 0;

            E_l_1 = 0;
            E_l_2_3 = 0;
            y_2_3 = 0;
            E_l_4 = 0;
            delta_4 = 0;
            E_l_5 = 0;
            E_l_6_7 = 0;
            E_l_8 = 0;

            cos_theta = 0.0; //0.99939082701909573000624344004393;
            sin_theta = 1.0; //0.034899496702500971645995181625333;

            //cos_theta := 0.98480775301220805936674302458952;
            //sin_theta := 0.17364817766693034885171662676931;

            Linear = false;

            SetAdditionalParameters();

            //AnnularMistTransitions(k1, k2);
        }
// Далее реализация интерфейса IHydrodynamicModel

        public void SetLinerWall(double innerDiameter, double relativeRoughness)
        {
            D = innerDiameter;
            h_pr = relativeRoughness;
        }

        public void SetInclination(double angle)
        {
            theta = angle;
        }

        public void SetHeavyFluid(double density, double viscosity)
        {
            rho_l = density;
            mu_l = viscosity;
        }

        public void SetLightFluid(double density, double viscosity)
        {
            rho_g = density;
            mu_g = viscosity;
        }

        public void SetSuperficialVelocityies(double light, double heavy)
        {
            V_sg = light;
            V_sl = heavy;
        }

        public void SetInterphaseTension(double tension)
        {
            sigma = tension;
        }

        public double GetHoldupLight()
        {
            return 1 - E_l;
        }

        public double GetHoldupHeavy()
        {
            return E_l;
        }

        public double GetPressureGradient()
        {
            return dP_dx;
        }

        public double GetVelocityLight()
        {
            return V_sg;   
        }

        public double GetVelocityHeavy()
        {
            return V_sl;
        }

        public void SetParameters(double p1, double p2, double p3, double p4, double p5, double p6, double p7)//Задание параметров и свойств для расчетов
                                                                                                       //те же параметры, что и в конструкторе
        {
            D = p1;
            rho_l = p2;
            rho_g = p3;
            mu_l = p4;
            mu_g = p5;
            sigma = p6;
            h_pr = p7;

            A = 0.25 * Math.PI * D * D;
            S = Math.PI * D;
        }

        public void SetConstParameters(double p1, double p2, double p3)//Задание константных параметров
        /*
                                                                          p1 - точность расчетов для расслоенного и кольцевого режимов течения
                                                                          p2 - ускорение свободного падения
                                                                          p3 - константа, необходимая для расчета расслоенного режима}
        */
        { 
            eps = p1;
            g = p2;
            sd = p3;
        }

        public void SetAdditionalParameters()//Расчет параметров, определяемых по свойствам системы
        {
            double c1, c2, c3, N_b;

            V_m = V_sl + V_sg;

            E_l_s = 1.0 / (1.0 + Math.Pow(V_m / 8.66, 1.39));

            N_b = (mu_l * V_sg / sigma) * (mu_l * V_sg / sigma) * rho_g / rho_l;
            c1 = 0.735 * Math.Pow(N_b, 0.074) * Math.Pow(V_sg / V_sl, 0.2);
            FE = c1 / (1.0 + c1);

            c2 = FE * V_sl / (FE * V_sl + V_sg);
            rho_c = c2 * rho_l + (1.0 - c2) * rho_g;
            mu_c = c2 * mu_l + (1.0 - c2) * mu_g;

            c3 = V_sg / (FE * V_sl + V_sg);
            delta_t = 0.5 * (1.0 - Math.Sqrt(0.76 / c3));
        }

        public double funH(double y_x, double V_sg_x, double V_sl_x)//Уравнение сохранения импульса (расслоенный режим)
        {
            double E_lsf,
            h_l_d,
            V_l, V_g, V_i,
            D_g,
            Re_g, Re_sl,
            Fr,
            f_sl, f_g, f_l, f_i,
            tau_wg, tau_wl, tau_i,
            S_l_d, S_g_d, S_i_d,
            A_l_d, A_g_d, A_d, A_S;

            S_l_d = y_x;        // S_l_d = 2*Pi*S_l/S, безразмерный периметр жидкой фазы, умноженный на 2*Pi
            S_g_d = 2.0 * Math.PI - S_l_d;
            S_i_d = 2.0 * Math.Sin(0.5 * y_x);
            A_l_d = y_x - Math.Sin(y_x);
            A_g_d = 2.0 * Math.PI - A_l_d;
            A_d = 2.0 * Math.PI;     // 2*Pi*A/A
            A_S = 0.25 * D;     // A/S размерная величина

            E_lsf = A_l_d / A_d;

            h_l_d = 0.5 * (1.0 - Math.Cos(0.5 * y_x));   //h_l_d = h_l/D

            V_l = V_sl_x / E_lsf;
            V_g = V_sg_x / (1.0 - E_lsf);
            V_i = V_g - V_l;

            D_g = A_S * GetHydravlicDiameter(A_g_d, S_g_d + S_i_d);

            Re_g = GetReynold(rho_g, D_g, V_g, mu_g);
            Re_sl = GetReynold(rho_l, D, V_sl_x, mu_l);
            Fr = V_l / Math.Sqrt(g * h_l_d * D);

            f_sl = GetFrictionFactor(Re_sl, D, h_pr);
            f_g = GetFrictionFactor(Re_g, D, h_pr);
            f_l = 0.452 * Math.Pow(f_sl, 0.731);
            f_i = (4e-3 + 5e-7 * Re_sl) * Math.Pow(Fr, 1.335) * (rho_l * D * g / (rho_g * V_g * V_g));

            tau_wg = GetShearStress(f_g, rho_g, V_g);
            tau_wl = GetShearStress(f_l, rho_l, V_l);
            tau_i = GetShearStress(f_i, rho_g, V_i);

            return (tau_wg * S_g_d * A_l_d - tau_wl * S_l_d * A_g_d + tau_i * S_i_d * A_d
                 - (rho_l - rho_g) * g * sin_theta * A_S * A_l_d * A_g_d);
        }

        public double funHy(double y_x)//Уравнение сохранения импульса (расслоенный режим)
        {
            return funH(y_x, V_sg, V_sl);
        }

        public double funD(double x, double V_sg_x, double V_sl_x)//Уравнение сохранения импульса (кольцевой-дисперсный режим)
        {
            double E_f, A_S,
            V_f, V_c,
            D_f, D_c, S_i, D_c_d,
            Re_f, Re_c,
            f_c, f_f, f_i,
            tau_wl, tau_i;

            E_f = 4.0 * x * (1.0 - x);
            A_S = 0.25 * D;

            V_f = V_sl_x * (1.0 - FE) / E_f;
            V_c = (V_sg_x + V_sl_x * FE) / (1.0 - E_f);

            D_c_d = (1.0 - 2.0 * x);
            D_c = D * D_c_d;
            //S_i   := Pi*D_c;
            D_f = GetHydravlicDiameter(A * E_f, S /* +S_i*/);

            Re_f = GetReynold(rho_l, D_f, V_f, mu_l);
            Re_c = GetReynold(rho_c, D_c, V_c, mu_c);
            if (Re_f < 1) { Re_f = 1.0; }
            if (Re_c < 1) { Re_c = 1.0; }

            f_c = GetFrictionFactor(Re_c, D, h_pr);
            f_f = GetFrictionFactor(Re_f, D, h_pr);
            f_i = f_c * 0.24 * Math.Pow(sigma / (rho_c * V_c * V_c * D_c), 0.085) * Math.Pow(Re_f, 0.305);

            tau_wl = GetShearStress(f_f, rho_l, V_f);
            tau_i = GetShearStress(f_i, rho_c, V_c - V_f);

            //из уравнения сохр-я импульса для пленки, умнож-ого на A_c, отнимается
            //уравнение сохр-я импульса для ядра, умнож-ое на A_f,
            //далее полученное уравнение умножается на 1/(A*S)
            return (-tau_wl * (1.0 - E_f) + tau_i * D_c_d
                     - A_S * E_f * (1.0 - E_f) * (rho_l - rho_c) * g * sin_theta);
        }

        public double funDx(double x)//Уравнение сохранения импульса (кольцевой-дисперсный режим)
        {
            return funD(x, V_sg, V_sl);
        }

        public double funDmin(double x)//Уравнение для расчета критической толщины пленки
        {
            double V_f,
            D_f,
            Re_f,
            f_f,
            E_f_x;

            E_f_x = 4.0 * x * (1.0 - x);

            V_f = V_sl * (1.0 - FE) / E_f_x;

            D_f = D * E_f_x;
            Re_f = GetReynold(rho_l, D_f, V_f, mu_l);
            // if Re_f<15 then Re_f := 15.0;

            f_f = GetFrictionFactor(Re_f, D, h_pr);

            //Result := E_f_x*sqr(E_f_x)*(1.0-1.5*E_f_x)/(2.0-1.5*E_f_x)
            //           - 2.0*f_f*rho_l*sqr(V_sl)*sqr(1.0-FE)/((rho_l-rho_c)*g*D*sin_theta);

            return (-f_f * rho_l * (V_sl * (1.0 - FE)) * (V_sl * (1.0 - FE)) * (4.0 - 3.0 * E_f_x) / (E_f_x * E_f_x * E_f_x)
                     + (rho_l - rho_c) * D * g * sin_theta * (1.0 - 1.5 * E_f_x));
        }

        public void DispersedBubble(ref Boolean Key_db, ref Boolean Key_theta, ref double HoldUp)
        {
            double C_g,
            Re_ml,
            C_o,
            V_b1, V_g_db,
            E_l_b;

            C_g = V_sg / V_m;

            Re_ml = rho_l * V_m * D / mu_l;                              //Co
            C_o = (1.64 + 0.12 * sin_theta) * Math.Pow(Re_ml, -0.031);

            V_b1 = 1.53 * Math.Pow((rho_l - rho_g) * g * sigma / (rho_l * rho_l), 0.25) * sin_theta;  //Vgdb
            V_g_db = C_o * V_m + V_b1;

            if (V_g_db > 0) { E_l_b = 1.0 - V_sg / V_g_db; } else { E_l_b = 1.0 - V_sg / (C_o * V_m); }
            if (E_l_b > 1) { E_l_b = V_sl / V_m; }

            Key_db = false;
            if ((E_l_s < 0.48) && (C_g <= 0.52)) { Key_db = true; }

            Key_theta = false;
            if (theta <= 0.087) { Key_theta = true; }

            HoldUp = E_l_b;
        }

        public void Stratified(ref Boolean Key_st, ref Boolean Key_st_sm, ref double HoldUp, ref double Gamma)
        {
            double y, h_l_d,
            E_l_sf,
            V_l, V_g,
            Re_sl,
            f_sl, f_l,
            A_g,
            dAldhl,
            cos_theta_s,
            V1_l, V1_g,
            V2_g, Fr_l,
            temp;

            y = SolveAlgEqN(eps, 2.0 * Math.PI - eps, eps, funHy, 100);

            h_l_d = 0.5 * (1.0 - Math.Cos(0.5 * y));
            E_l_sf = (y - Math.Sin(y)) / (2.0 * Math.PI);

            /*
            h_l_d:= SolveAlgEqN(eps, 1.0 - eps, eps, funHy, 100);
            temp:= 2.0 * h_l_d - 1.0;
            E_l_sf:= (pi - arccos(temp) + temp * sqrt(1.0 - sqr(temp))) / pi;
            */

            V_l = V_sl / E_l_sf;
            V_g = V_sg / (1.0 - E_l_sf);

            Re_sl = GetReynold(rho_l, D, V_l, mu_l);
            f_sl = GetFrictionFactor(Re_sl, D, h_pr);
            f_l = 0.452 * Math.Pow(f_sl, 0.731);

            A_g = A * (1.0 - E_l_sf);
            dAldhl = 2.0 * A * (1.0 - Math.Cos(y)) / (Math.PI * D * Math.Sin(0.5 * y));
            //dAldhl   := D*sqrt(1.0 - sqr(temp));
            cos_theta_s = cos_theta;
            if (cos_theta_s < 0.02) { cos_theta_s = 0.02; }

            V1_l = Math.Sqrt(g * D * (1.0 - h_l_d) * cos_theta / f_l);
            V1_g = (1.0 - h_l_d) * Math.Sqrt((rho_l - rho_g) * g * A_g * cos_theta_s / (rho_g * dAldhl));

            V2_g = Math.Sqrt(4.0 * mu_l * (rho_l - rho_g) * g * cos_theta / (sd * rho_l * rho_g * V_l));
            Fr_l = V_l / Math.Sqrt(g * h_l_d * D);

            Key_st = false;
            if ((V_l <= V1_l) && (V_g <= V1_g)) { Key_st = true; }

            Key_st_sm = false;
            if ((V_g <= V2_g) && (Fr_l <= 1.4)) { Key_st_sm = true; }

            Gamma = y;
            //Gamma  := h_l_d;
            HoldUp = E_l_sf;
        }

        public void AnnularMist(ref Boolean Key_am, ref Boolean Key_f, ref double HoldUp, ref double Thickness)
        {
            double delta_d,
            delta_min_d,
            E_l_am;

            delta_d = SolveAlgEqN(eps, 0.5 - eps, eps, funDx, 100);


            /*
                SetLength(fun, 1001);
                for i := 0 to 1000 do
                        SetLength(fun[i], 2);



                fun[0, 0] := eps;
                fun[1000, 0] := 0.5 - eps;
                fun[0, 1] := funDx(eps);
                fun[1000, 1] := funDx(0.5 - eps);
                for i := 1 to 999 do
                        begin
                          fun[i, 0] := 0.5 * i / 1000;
                fun[i, 1] := funDx(0.5 * i / 1000);
                end;

                SaveMatrixData(fun, 'Result.txt');
            */

            E_l_am = 1.0 - (1.0 - 2.0 * delta_d) * (1.0 - 2.0 * delta_d) * V_sg / (V_sg + FE * V_sl);

            Thickness = delta_d;
            HoldUp = E_l_am;

            if (theta > 0)
            {
                delta_min_d = SolveAlgEqN(eps, 0.5 - eps, eps, funDmin, 100);
                if (double.IsNaN(delta_min_d)) { delta_min_d = 1.0; }
            }
            else
            { delta_min_d = 1.0; }

            Key_am = false;
            if ((E_l_am <= 0.24) /* (delta_d <= delta_t)*/ && (delta_d < delta_min_d)) { Key_am = true; }

            Key_f = false;
            if (E_l_s < 0.48) { Key_f = true; }
        }

        public void BubbleAndIntermitten(ref Boolean Key_b, ref Boolean Key_s_eb, ref Boolean Key_s,
                                  ref double HoldUp_1, ref double HoldUp_2)
        {
            double V_b,
            Re_ml,
            C_o, C_l, gamma,
            Bo, V_d_h_inf, beta, V_d_nu_inf, V_d_inf,
            Re_inf, f_m, V_d,
            V_t,
            V_b1, V_g_db,
            D1, D_b,
            cos1_theta,
            E_l_slug, E_g, C_g;

            D1 = 19.0 * Math.Sqrt((rho_l - rho_g) * sigma / (rho_l * rho_l * g));

            V_b = 1.41 * Math.Pow((rho_l - rho_g) * g * sigma / (rho_l * rho_l), 0.25) * sin_theta;
            C_l = 0.8;
            gamma = 1.3;
            D_b = 7e-3;
            cos1_theta = 0.53003 * V_b * V_b * (C_l * gamma * gamma / (g * D_b));

            Re_ml = GetReynold(rho_l, D, V_m, mu_l);                                     //Co для пробкового
            C_o = (1.64 + 0.12 * sin_theta) * Math.Pow(Re_ml, -0.031);

            Bo = (rho_l - rho_g) * g * D * D / sigma;                                    //Vdinf
            V_d_h_inf = (0.54 - 1.76 / Math.Pow(Bo, 0.56)) * Math.Sqrt(g * D * (rho_l - rho_g) / rho_l);
            beta = Bo * Math.Exp(3.278 - 1.424 * Math.Log(Bo));
            V_d_nu_inf = 0.345 * (1.0 - Math.Exp(-beta)) * Math.Sqrt(g * D * (rho_l - rho_g) / rho_l);
            V_d_inf = V_d_h_inf * cos_theta + V_d_nu_inf * sin_theta;

            Re_inf = rho_l * V_d_inf * D / (2.0 * mu_l);                                         //Vd
            f_m = 0.316 * Math.Sqrt(Re_inf);
            if (f_m > 1) { f_m = 1.0; }
            V_d = f_m * V_d_inf;

            V_t = C_o * V_m + V_d;                                                         //Vt для продолговато-пузырькового

            V_b1 = 1.53 * Math.Pow((rho_l - rho_g) * g * sigma / (rho_l * rho_l), 0.25) * sin_theta;    //Vgdb
            V_g_db = C_o * V_m + V_b1;
            if (V_g_db < 0) { V_g_db = 0; }

            E_l_slug = (E_l_s * V_t + V_g_db * (1.0 - E_l_s) - V_sg) / V_t;                          //Elslug
            if (E_l_slug > 1) { E_l_slug = V_sl / V_m; }

            C_o = 1.2;                                                                   //Co для пузырькового
            V_t = C_o * V_m + V_b;                                                         //Vt для пузырькового
            E_g = V_sg / V_t;                                                              //Eg
            C_g = V_sg / V_m;
            if (E_g > C_g) { E_g = C_g; }

            Key_b = false;
            if ((cos_theta <= cos1_theta) && (D > D1) && (E_l_slug > 0.75)) { Key_b = true; }

            Key_s_eb = false;
            if (E_l_slug > 0.24) { Key_s_eb = true; }

            Key_s = false;
            if (E_l_s <= 0.9) { Key_s = true; }

            HoldUp_1 = 1.0 - E_g;
            HoldUp_2 = E_l_slug;
        }

        public void DispersedBubbleTransitions(ref Boolean Key_db, ref Boolean Key_theta)//Расчет фазового содержания и переходов (дисперсный режим)
        {
            DispersedBubble(ref Key_db, ref Key_theta, ref E_l_1);
        }

        public void StratifiedTransitions(ref Boolean Key_st, ref Boolean Key_st_sm)//Расчет фазового содержания и переходов (расслоенный режим)
        {
            Stratified(ref Key_st, ref Key_st_sm, ref E_l_2_3, ref y_2_3);
        }

        public void AnnularMistTransitions(ref Boolean Key_am, ref Boolean Key_f)//Расчет фазового содержания и переходов (кольцевой-дисперсный режим)
        {
            AnnularMist(ref Key_am, ref Key_f, ref E_l_4, ref delta_4);
        }

        public void BubbleAndIntermittenTransitions(ref Boolean Key_b, ref Boolean Key_s_eb, ref Boolean Key_s)//Расчет фазового содержания и переходов (пузырьковый, снарядный и пробковый режимы)
        {
            BubbleAndIntermitten(ref Key_b, ref Key_s_eb, ref Key_s, ref E_l_5, ref E_l_6_7);
        }

        public void Calculate()//Расчет режима течения по заданным расходам
        {
            Boolean Key_db = true, Key_theta = true, Key_st = true, Key_st_sm = true, Key_am = true, Key_f = true, Key_b = true, Key_s_eb = true, Key_s = true;

            /*
              DispersedBubble  = 1
              StratifiedSmooth = 2
              StratifiedWavy   = 3
              AnnularMist      = 4
              Bubble           = 5
              Slug             = 6
              ElongatedBubble  = 7
              Froth            = 8
            */
            SetAdditionalParameters();

            DispersedBubbleTransitions(ref Key_db, ref Key_theta);

            if (Key_db)
            {
                FlowRegime = 1;
                return;
            }

            if (Key_theta)
            {
                StratifiedTransitions(ref Key_st, ref Key_st_sm);
                if (Key_st)
                {
                    if (Key_st_sm)
                    {
                        FlowRegime = 2;
                        return;
                    }
                    else
                    {
                        FlowRegime = 3;
                        return;
                    }
                }
            }
                    AnnularMistTransitions(ref Key_am, ref Key_f);
                    if (Key_am)
                    {
                        FlowRegime = 4;
                        return;
                    }
                    else
             if (Key_f)
                    {
                        FlowRegime = 8;
                        return;
                    }
                    else
                    {
                        BubbleAndIntermittenTransitions(ref Key_b, ref Key_s_eb, ref Key_s);
                        if (Key_b)
                        {
                            FlowRegime = 5;
                            return;
                        }
                        else
                        {
                            if (Key_s_eb)
                            {
                                if (Key_s)
                                {
                                    FlowRegime = 6;
                                    return;
                                }
                                else
                                {
                                    FlowRegime = 7;
                                    return;
                                }
                            }
                            else
                            {
                                FlowRegime = 8;
                                return;
                            }
                        }
                    }
                }

        public void CalcLiquidVolumeFraction()//Расчет фазового содержания жидкости
        {
            Calculate();
            if ((FlowRegime == 1) || (FlowRegime == 8)) { E_l = E_l_1; }     //!!!!!!!!!!!!!
            if ((FlowRegime == 2) || (FlowRegime == 3)) { E_l = E_l_2_3; }
            if (FlowRegime == 4) { E_l = E_l_4; }
            if (FlowRegime == 5) { E_l = E_l_5; }
            if ((FlowRegime == 6) || (FlowRegime == 7)) { E_l = E_l_6_7; }
        }

        public void CalcPressureGradient()//Расчет градиента давления
        {
            double rho_m, mu_m,
            Re_m, Re_ml, Re_f, Re_c, Re_sl, Fr,
            f_c, f_i, f_f, f_m1, f_m2, f_sl, f_l,
            S_l, S_i, S_c,
            D_f, D_c, D_c_d,
            V_l, V_g, V_f, V_c, V_i,
            tau_wl, tau_i,
            h_l, delta_i, eta,
            E_f, Sc_Ac;

            Calculate();
            CalcLiquidVolumeFraction();

            if ((FlowRegime == 1) || (FlowRegime == 5))
            {
                rho_m = rho_g * (1.0 - E_l) + rho_l * E_l;
                mu_m = mu_g * (1.0 - E_l) + mu_l * E_l;

                Re_m = GetReynold(rho_m, D, V_m, mu_m);
                if (Re_m < 15) { Re_m = 15.0; }
                f_m1 = GetFrictionFactor(Re_m, D, h_pr);

                dP_dx = 2.0 * f_m1 * V_m * V_m * rho_m / D + rho_m * g * sin_theta;
            }

            if ((FlowRegime == 2) || (FlowRegime == 3))
            {
                h_l = 0.5 * D * (1.0 - Math.Cos(0.5 * y_2_3));

                S_l = 0.5 * y_2_3 * D;
                S_i = Math.Sin(0.5 * y_2_3) * D;

                V_l = V_sl / E_l;
                V_g = V_sg / (1.0 - E_l);
                V_i = V_g - V_l;

                Re_sl = GetReynold(rho_l, D, V_sl, mu_l);
                Fr = V_l / Math.Sqrt(g * h_l);

                f_sl = GetFrictionFactor(Re_sl, D, h_pr);
                f_l = 0.452 * Math.Pow(f_sl, 0.731);
                f_i = (4e-3 + 5e-7 * Re_sl) * Math.Pow(Fr, 1.335) * (rho_l * D * g / (rho_g * V_g * V_g));

                tau_wl = GetShearStress(f_l, rho_l, V_l);
                tau_i = GetShearStress(f_i, rho_g, V_i);

                dP_dx = (tau_wl * S_l - tau_i * S_i) / (E_l * A) + rho_l * g * sin_theta;
            }

            if (FlowRegime == 4)
            {
                E_f = 4.0 * delta_4 * (1.0 - delta_4);

                V_f = V_sl * (1.0 - FE) / E_f;
                V_c = V_sg / (1.0 - E_f);

                D_c_d = (1.0 - 2.0 * delta_4);
                D_c = D * D_c_d;
                S_c = Math.PI * D_c;
                D_f = GetHydravlicDiameter(A * E_f, S + S_c);
                Sc_Ac = 4.0 / D_c;

                Re_f = GetReynold(rho_l, D_f, V_f, mu_l);
                Re_c = GetReynold(rho_c, D_c, V_c, mu_c);
                if (Re_f < 15) { Re_f = 15.0; }
                if (Re_c < 15) { Re_c = 15.0; }

                f_c = GetFrictionFactor(Re_c, D, h_pr);
                f_i = f_c * 0.24 * Math.Pow(sigma / (rho_c * V_c * V_c * D_c), 0.085) * Math.Pow(Re_f, 0.305);

                tau_i = GetShearStress(f_i, rho_c, V_c - V_f);

                dP_dx = tau_i * Sc_Ac + rho_c * g * sin_theta;
            }


            if ((FlowRegime == 6) || (FlowRegime == 7))
            {
                rho_m = rho_g * (1.0 - E_l) + rho_l * E_l;
                mu_m = mu_g * (1.0 - E_l) + mu_l * E_l;

                delta_i = 0.5 * (1 - Math.Sqrt((1.0 - E_l) * (FE * V_sl + V_sg) / V_sg));

                eta = Math.Pow(V_sl / V_m, 0.75 - E_l);

                if (eta > 1) { eta = 1.0; }


                Re_ml = GetReynold(rho_l, D, V_m, mu_l);
                f_m1 = GetFrictionFactor(Re_ml, D, h_pr);

                if (delta_i >= 1e-4)
                {
                    E_f = 4.0 * delta_i * (1.0 - delta_i);
                    V_f = V_sl * (1 - FE) / E_f;
                    D_c_d = (1.0 - 2.0 * delta_4);
                    D_c = D * D_c_d;
                    S_c = Math.PI * D_c;
                    D_f = GetHydravlicDiameter(A * E_f, S + S_c);

                    Re_f = GetReynold(rho_l, D_f, V_f, mu_l);
                    f_f = GetFrictionFactor(Re_f, D, h_pr);
                    tau_wl = GetShearStress(f_f, rho_l, V_f);

                    dP_dx = rho_m * g * sin_theta + eta * (2 * f_m1 * V_m * V_m * rho_m / D) + (1.0 - eta) * (4.0 * tau_wl / D);
                }
                else
                {
                    Re_m = GetReynold(rho_m, D, V_m, mu_m);
                    f_m2 = GetFrictionFactor(Re_m, D, h_pr);

                    dP_dx = rho_m * g * sin_theta + eta * (2 * f_m1 * V_m * V_m * rho_m / D) + (1.0 - eta) * (2.0 * f_m2 * V_m * V_m * rho_m / D);
                }
            }
            if (FlowRegime == 8)
            {
                dP_dx = 0;

            }
        }

        public double funEta(double H_l)//Функция, показывающая степень влияния каждой фазы на показания расходомера (!!!)
        {
            double Eta;

            Eta = 0.0;

            if (H_l >= 0.7) { Eta = 1.0; }
            if (H_l <= 0.3) { Eta = 0.0; }
            if (Linear)
            {
                if ((H_l < 0.7) && (H_l > 0.3)) { Eta = 2.5 * (H_l - 0.3); }
            } else {
                if ((H_l < 0.7) && (H_l > 0.3)) { Eta = 1.015 * (1 - Math.Exp(-10.0 * (H_l - 0.3))); }
            }

        return Eta;
        }

        public double funEl_Vsg(double V_x, double V, double H_l)
        {
            double Eta;

        Eta = funEta(H_l);

        V_sg = V_x;
        V_sl = (V - (1 - Eta) * V_x) / Eta;

            CalcLiquidVolumeFraction();

        return E_l - H_l;
        }

        public double funEl_Vsl(double V_x, double V, double H_l)
        {
            double Eta;

            Eta = funEta(H_l);
            V_sl = V_x;
            V_sg = (V - Eta * V_x) / (1 - Eta);

            CalcLiquidVolumeFraction();

        return E_l - H_l;
        }

        public double CalcVSG(double a, double b, double V, double H_l)
        {
            double c;

        c = 0;

            do
            {
                c = (b + a) / 2;
                if (funEl_Vsg(a, V, H_l) * funEl_Vsg(c, V, H_l) < 0) { b = c; } else { a = c; }
            } while (Math.Abs(b - a) > eps);

        return c;
        }

        public double CalcVSL(double a, double b, double V, double H_l)
        {
            double c;

            c = 0;

            do
            {
                c = (b + a) / 2;
                if (funEl_Vsl(a, V, H_l) * funEl_Vsl(c, V, H_l) < 0) { b = c; } else { a = c; }
            } while (Math.Abs(b - a) > eps);

            return c;
        }

        public void GetVelocities(double Velocity, double LiquidHoldUp, double Angle, ref double SLVelocity, ref double SGVelocity)//Расчет приведенных скоростей (расходов) по заданному фазовому содержанию и видимой скорости
        {
            double V_sl_calc, V_sg_calc, Eta, V;

        theta = Angle * Math.PI / 180;
        cos_theta = Math.Cos(theta);
        sin_theta = Math.Sin(theta);

        V_sl_calc = 0;
        V_sg_calc = 0;

            if (LiquidHoldUp >= 0.5)
            {
                Eta = funEta(LiquidHoldUp);
                if (Eta == 1) { V = 1000; }
                else { V = Velocity / (1 - Eta);
                    if (V > 1000) { V = 1000; }
                }
                V_sg_calc = CalcVSG(eps, V - eps, Velocity, LiquidHoldUp);
                V_sl_calc = (Velocity - (1 - Eta) * V_sg_calc) / Eta;
            }

            if (LiquidHoldUp < 0.5)
            {
                Eta = funEta(LiquidHoldUp);
                if (Eta == 0) { V = 1000; }
                else { V = Velocity / Eta;
                    if (V > 1000) { V = 1000; }
                }
                V_sl_calc = CalcVSL(eps, V - eps, Velocity, LiquidHoldUp);
                V_sg_calc = (Velocity - Eta * V_sl_calc) / (1 - Eta);
            }

        SLVelocity = V_sl_calc;
        SGVelocity = V_sg_calc;
        }

        public byte GetFlowRegime(double Velocity_SL, double Velocity_SG, double Angle)//Расчет режима течения по заданным приведенным скоростям (расходам)
        {
            theta = Angle * Math.PI / 180;
            cos_theta = Math.Cos(theta);
            sin_theta = Math.Sin(theta);

            V_sl = Velocity_SL;
            V_sg = Velocity_SG;

                Calculate();

            return FlowRegime;
        }

        public double GetLiquidHoldup(double Velocity_SL, double Velocity_SG, double Angle)//Расчет фазового содержания жидкости по заданным приведенным скоростям (расходам)
        {
            theta = Angle * Math.PI / 180;
            cos_theta = Math.Cos(theta);
            sin_theta = Math.Sin(theta);

            V_sl = Velocity_SL;
            V_sg = Velocity_SG;

            CalcLiquidVolumeFraction();

            return E_l;
        }

        public double GetPressureGradient(double Velocity_SL, double Velocity_SG, double Angle)//Расчет градиента давления по заданным приведенным скоростям (расходам)
        {
            theta = Angle * Math.PI / 180;
            cos_theta = Math.Cos(theta);
            sin_theta = Math.Sin(theta);

            V_sl = Velocity_SL;
            V_sg = Velocity_SG;

            CalcPressureGradient();

            return dP_dx;
        }

        public double GetHydravlicDiameter(double Area, double Perimeter)
        {
            if ((Area<=0) || (Perimeter<=0))
            {
                throw new System.ArgumentException("Отрицательные или нулевые величины входных параметров (GetShearStress)");
            }
            return 4.0 * Area / Perimeter;
        }

        public double GetReynold(double Density, double Length, double Velocity, double Viscosity)
        {
            if ((Density <= 0) || (Length<=0) || (Viscosity<=0))
            {
                throw new System.ArgumentException("Отрицательные или нулевые величины входных параметров (GetReynold)");   
            }
            return Density * Length * Math.Abs(Velocity) / Viscosity;
        }

        public double GetFrictionFactor(double Reynold, double PipeDiameter, double PipeRoughness)
        {
            double x, e, re_1, f, ft, fl;

            if ((Reynold <= 0) || (PipeDiameter <= 0) || (PipeRoughness < 0))
            {
                throw new System.ArgumentException("Отрицательные или нулевые величины входных параметров (GetFrictionFactor)");
            }

            if (PipeRoughness >= PipeDiameter)
            {
                throw new System.ArgumentException("Величина абсолютной неровности трубы не может быть больше диаметра трубы (GetFrictionFactor)");
            }

            try
            {
                re_1 = 1.0 / Reynold;
                if (Reynold < 20)
                { f = 16.0 * re_1; }
                else
                {
                    e = PipeRoughness / PipeDiameter;
                    x = -4.0 * Math.Log10(0.2698 * e - 5.0452 * re_1 * Math.Log10(0.3539 * Math.Pow(e, 1.1098) + 5.8506 * Math.Pow(re_1, 0.8981)));
                    ft = 1.0 / x / x;
                    if (Reynold < 2000)
                    {
                        fl = 16.0 * re_1;
                        if (ft > fl) { f = ft; } else { f = fl; }

                    }
                    else
                        f = ft;
                }
                return f;
            }

            catch
            {
                throw new Exception("Ошибка при расчете коэффициента трения (GetFrictionFactor)");
            }
        }

        public double GetShearStress(double FrictionFactor, double Density, double Velocity)
        {
            if ((FrictionFactor < 0) || (Density <= 0))
            {
                throw new ArgumentException("Отрицательные или нулевые величины входных параметров (GetShearStress)");
            }
            return 0.5 * FrictionFactor * Density * Velocity * Math.Abs(Velocity);
        }

        public double SolveAlgEq(double LBound, double RBound, double Err, Func<double, double> F, double FuncRes )
        {
            double MaxIter = 1e8;

            int i;
            double X_C, X_C_pr, X_L, X_R, Y_C, Y_L, Y_R;

            X_L = LBound;
            X_R = RBound;
            X_C = LBound;
            Y_L = F(X_L) - FuncRes;
            Y_R = F(X_R) - FuncRes;

            if (Y_L * Y_R > 0) { throw new ArgumentException("На отрезке нет корней или их количество больше одного (четное) (SolveAlgEq)"); }

            i = 0;

            do
            {
                X_C_pr = X_C;
                X_C = X_L - Y_L * (X_R - X_L) / (Y_R - Y_L);
                Y_C = F(X_C) - FuncRes;
                if (Math.Abs(X_C - X_C_pr) < Err)
                {
                    return X_C;
                }
                if (Y_L * Y_C < 0)
                {
                    X_R = X_C;
                    Y_R = Y_C;
                }
                else
                {
                    X_L = X_C;
                    Y_L = Y_C;
                }
                i++;
                if (i == MaxIter) { throw new Exception("Решение уравнения не сходится (SolveAlgEq)"); }
            }
            while (Math.Abs(X_C - X_C_pr) > Err);
            return Y_C;
        }

        public double SolveAlgEqN(double LBound, double RBound, double Err, Func<double, double> F, double FuncRes, int N)
        {
            double dx, a1, a, b;
            int i;

            if (N == 0)
            {
                throw new Exception("N не может равняться 0 (SolveAlgEq)");
            }
            dx = (RBound - LBound) / N;
            a1 = LBound;

            for (i = 0; i < (N - 1); i++)
            {
                a = a1 + i * dx;
                b = a + dx;
                if (F(a) * F(b) < 0)
                {
                    return SolveAlgEq(a, b, Err, F, FuncRes);
                } //if
            }  //for
        return double.NaN;
        }

        public double SolveAlgEqN(double LBound, double RBound, double Err, Func<double, double> F, int N)
        {
            return SolveAlgEqN(LBound, RBound, Err, F, 0, N);
        }

    }
}
