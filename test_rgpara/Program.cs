using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data;
using System.Runtime.CompilerServices;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Globalization;
//using MathNet.Numerics.LinearAlgebra.Generic;

namespace test_rgpara
{
    class Program
    {
        static void Main(string[] args)
        {
            //test_rgpara();
            test_rgpara2();
        }
        /// <summary>
        /// 计算全息特征参数(调用MathNet.Numerics库，对向量进行操作)
        /// </summary>
        public static void rgpara(Vector<double> VectorR, out double Long_a, out double Short_b, out double E, out double Fi)
        {
            var M = Matrix<double>.Build;
            var V = Vector<double>.Build;
            var formatProvider = (CultureInfo)CultureInfo.InvariantCulture.Clone();
            formatProvider.TextInfo.ListSeparator = " ";

            //目的：椭圆参数转换
            //输入变量：r=[sx cx sy cy];
            var SX = VectorR[0];//  sx - 椭圆x分量的正弦系数;
            var CX = VectorR[1];//  cx －椭圆x分量的余弦系数;
            var SY = VectorR[2];//  sy - 椭圆y分量的正弦系数;
            var CY = VectorR[3];//  cy - 椭圆y分量的余弦系数;

            //输出变量：
            //a,b 为椭圆半长轴与半短轴;
            //fi: 椭圆长轴的倾角;
            //e:椭圆偏心率;

            double R2a = Math.PI / 180; //const=pi/180;
            //t=1:361;
            //x=r(1).*sin((t-1).*const)+r(2).*cos((t-1).*const);
            //y=r(3).*sin((t-1).*const)+r(4).*cos((t-1).*const);
            Vector<double> VectorCircle= V.Dense(361, i => R2a * i);
            // Format vector output to console
            //Console.WriteLine(@"VectorCircle = ");
            //Console.WriteLine(VectorCircle.ToString("#0.00\t", formatProvider));
            //Console.WriteLine();

            Vector<double> vectorCircleSin = V.Dense(361, i => Math.Sin(VectorCircle[i]));
            Vector<double> vectorCircleCos = V.Dense(361, i => Math.Cos(VectorCircle[i]));

            var VectorX = SX * vectorCircleSin + CX * vectorCircleCos;
            var VectorY = SY * vectorCircleSin + CY * vectorCircleCos;
            //Console.WriteLine(@"VectorX = ");
            //Console.WriteLine(VectorX.ToString("#0.00\t", formatProvider));
            //Console.WriteLine();
            //Console.WriteLine(@"VectorY = ");
            //Console.WriteLine(VectorY.ToString("#0.00\t", formatProvider));
            //Console.WriteLine();

            var VectorTanXY = VectorX / VectorY;//tan_xy=y./x;
            //Console.WriteLine(@"VectorTanXY= ");
            //Console.WriteLine(VectorTanXY.ToString("#0.00\t", formatProvider));
            //Console.WriteLine();

            //double p = Math.Pow(sx, 2.0) + Math.Pow(sy, 2.0) + Math.Pow(cx, 2.0) + Math.Pow(cy, 2.0);
            double P = Math.Pow(VectorR.L2Norm(), 2.0); //p = r(1).^2 + r(2).^2 + r(3).^2 + r(4).^2;
            double Q = 2 * CX * SY + 2 * SX * CY;   //q = 2.*r(2).*r(3) - 2.*r(1).*r(4);

            //输出1 2
            Long_a  =          (Math.Sqrt(P + Q) + Math.Sqrt(P - Q)) / 2;//long_a=(sqrt(p+q)+sqrt(p-q))/2;
            Short_b = Math.Abs((Math.Sqrt(P + Q) - Math.Sqrt(P - Q)) / 2);//short_b=abs((sqrt(p+q)-sqrt(p-q))/2);

            //输出3
            E = Math.Sqrt( Math.Pow(Long_a, 2.0) - Math.Pow(Short_b, 2.0) ) / Long_a;//e = sqrt(long_a^2 -short_b^2)/long_a;
            //SC=(2*r(2)*r(4)+2*r(1)*r(3))./(r(2).^2+r(1).^2-r(4).^2-r(3).^2);
            double SC = (2 * VectorR[1] * VectorR[3] + 2 * VectorR[0] * VectorR[1]) / (Math.Pow(SX, 2.0) + Math.Pow(CX, 2.0) - Math.Pow(SY, 2.0) - Math.Pow(CY, 2.0));
            double Angle_d = 0.5 * Math.Atan(SC) / R2a;//angle_d=0.5*atan(SC)./const;
            double Tan_angle = Math.Tan(Angle_d * R2a);//tan_angle=tan(angle_d.*const);
            Vector<double> Cmp = V.Dense(361, i => Math.Abs( VectorTanXY[i] - Tan_angle) );//cmp=abs(tan_xy-tan_angle);

            int Mino = Cmp.MinimumIndex();//[mincmp,mino]=min(cmp);
            double Axis_cmp = Math.Pow(VectorX[Mino], 2.0) + Math.Pow(VectorY[Mino], 2.0);//axis_cmp=x(mino).^2+y(mino).^2;
            double ACmp = Math.Abs(Math.Pow(Long_a, 2.0) - Axis_cmp);//acmp=abs(long_a.^2-axis_cmp);
            double BCmp = Math.Abs(Math.Pow(Short_b, 2.0) - Axis_cmp);//bcmp=abs(short_b.^2-axis_cmp);
            //if acmp>bcmp
            //    angle_d=angle_d+90;
            //end  
            if(ACmp > BCmp)
            {
                Angle_d += 90;
            }
            //if angle_d<0
            //    angle_d=angle_d+180;
            //end
            if (Angle_d < 0)
            {
                Angle_d += 180;
            }
            //if angle_d>=180
            //    angle_d=angle_d-180;
            //end
            if (Angle_d >= 180)
            {
                Angle_d -= 180;
            }
            //输出4
            Fi = Angle_d;//fi=angle_d;
            //IPP=[r(2) r(4)];
        }
        static void test_rgpara()
        {
            //test函数rgpara
            var M = Matrix<double>.Build;
            var V = Vector<double>.Build;
            var formatProvider = (CultureInfo)CultureInfo.InvariantCulture.Clone();
            formatProvider.TextInfo.ListSeparator = " ";
            var vectorR = V.Dense(new[] { -3.1829231e+01, 5.8633108e+00, -3.3810352e+00, -3.8530371e+01 });
            //Console.WriteLine(@"vectorR = ");
            //Console.WriteLine(vectorR.ToString("#0.00\t", formatProvider));
            //Console.WriteLine();
            //声明输出
            double Long_a;
            double Short_b;
            double E;
            double Fi;
            rgpara(vectorR, out Long_a, out Short_b, out E, out Fi);
            Console.WriteLine("计算结果如下: ");
            Console.WriteLine("Long_a = {0}", Long_a);
            Console.WriteLine("Short_b = {0}", Short_b);
            Console.WriteLine("E = {0}", E);
            Console.WriteLine("Fi = {0}", Fi);
            Console.ReadKey();
        }
        /// <summary>
        /// 计算全息特征参数(使用循环操作数组)
        /// </summary>
        public static void rgpara2(double[] ArrayR, out double Long_a, out double Short_b, out double E, out double Fi)
        {
            //目的：椭圆参数转换
            //输入变量：r=[sx cx sy cy];
            var Sx = ArrayR[0];//  sx - 椭圆x分量的正弦系数;
            var Cx = ArrayR[1];//  cx －椭圆x分量的余弦系数;
            var Sy = ArrayR[2];//  sy - 椭圆y分量的正弦系数;
            var Cy = ArrayR[3];//  cy - 椭圆y分量的余弦系数;

            //输出变量：
            //a,b 为椭圆半长轴与半短轴;
            //fi: 椭圆长轴的倾角;
            //e:椭圆偏心率;

            double R2a = Math.PI / 180; //const=pi/180;

            //t=1:361;
            //x=r(1).*sin((t-1).*const)+r(2).*cos((t-1).*const);
            //y=r(3).*sin((t-1).*const)+r(4).*cos((t-1).*const);
            double[] ArrayCircle = new double[361];
            double[] ArrayCircleSin = new double[361];
            double[] ArrayCircleCos = new double[361];
            //中间变量
            double[] SxByArrayCircleSin = new double[361];//Sx * ArrayCircleSin
            double[] CxByArrayCircleCos = new double[361];//Cx * ArrayCircleCos
            double[] SyByArrayCircleSin = new double[361];//Sy * ArrayCircleSin
            double[] CyByArrayCircleCos = new double[361];//Cy * ArrayCircleCos

            double[] ArrayX = new double[361];
            double[] ArrayY = new double[361];
            double[] ArrayTanXY = new double[361];

            for (int i = 0; i < ArrayCircle.Length; i++)
            {
                ArrayCircle[i] = i;
                ArrayCircleSin[i] = Math.Sin(i);
                ArrayCircleCos[i] = Math.Cos(i);

                SxByArrayCircleSin[i] = Sx * ArrayCircleSin[i];
                CxByArrayCircleCos[i] = Cx * ArrayCircleCos[i];
                SyByArrayCircleSin[i] = Sy * ArrayCircleSin[i];
                CyByArrayCircleCos[i] = Cy * ArrayCircleCos[i];

                ArrayX[i] = SxByArrayCircleSin[i] + CxByArrayCircleCos[i];
                ArrayY[i] = SyByArrayCircleSin[i] + CyByArrayCircleCos[i];
                ArrayTanXY[i] = ArrayX[i] / ArrayY[i];
            }

            double P = Math.Pow(Sx, 2.0) + Math.Pow(Sy, 2.0) + Math.Pow(Cx, 2.0) + Math.Pow(Cy, 2.0);//p = r(1).^2 + r(2).^2 + r(3).^2 + r(4).^2;
            double Q = 2 * Cx * Sy + 2 * Sx * Cy;   //q = 2.*r(2).*r(3) - 2.*r(1).*r(4);

            //输出1 2
            Long_a = (Math.Sqrt(P + Q) + Math.Sqrt(P - Q)) / 2;//long_a=(sqrt(p+q)+sqrt(p-q))/2;
            Short_b = Math.Abs((Math.Sqrt(P + Q) - Math.Sqrt(P - Q)) / 2);//short_b=abs((sqrt(p+q)-sqrt(p-q))/2);

            //输出3
            E = Math.Sqrt(Math.Pow(Long_a, 2.0) - Math.Pow(Short_b, 2.0)) / Long_a;//e = sqrt(long_a^2 -short_b^2)/long_a;
            //SC=(2*r(2)*r(4)+2*r(1)*r(3))./(r(2).^2+r(1).^2-r(4).^2-r(3).^2);
            double SC = (2 * ArrayR[1] * ArrayR[3] + 2 * ArrayR[0] * ArrayR[1]) / (Math.Pow(Sx, 2.0) + Math.Pow(Cx, 2.0) - Math.Pow(Sy, 2.0) - Math.Pow(Cy, 2.0));
            double Angle_d = 0.5 * Math.Atan(SC) / R2a;//angle_d=0.5*atan(SC)./const;
            double Tan_angle = Math.Tan(Angle_d * R2a);//tan_angle=tan(angle_d.*const);

            //cmp = abs(tan_xy - tan_angle);
            double[] Cmp = new double[361];
            for (int i = 0; i < ArrayCircle.Length; i++)
            {
                //Cmp[i] = Math.Cos(ArrayTanXY[i]) - Tan_angle;
                Cmp[i] = Math.Abs( ArrayTanXY[i] - Tan_angle );
            }
            
            //[mincmp,mino]=min(cmp);
            double CmpMin = Cmp[0];
            int CmpMino = 0;//把假设的最大值索引赋值非index
            for (int i = 1; i < Cmp.Length; i++)
            {
                if (Cmp[i] < CmpMin)
                {
                    CmpMin = Cmp[i];
                    CmpMino = i;//把较小值的索引赋值非index
                }
            }
            double Axis_cmp = Math.Pow(ArrayX[CmpMino], 2.0) + Math.Pow(ArrayY[CmpMino], 2.0);//axis_cmp=x(mino).^2+y(mino).^2;
            double ACmp = Math.Abs(Math.Pow(Long_a, 2.0) - Axis_cmp);//acmp=abs(long_a.^2-axis_cmp);
            double BCmp = Math.Abs(Math.Pow(Short_b, 2.0) - Axis_cmp);//bcmp=abs(short_b.^2-axis_cmp);

            //if acmp>bcmp
            //    angle_d=angle_d+90;
            //end  
            if (ACmp > BCmp)
            {
                Angle_d += 90;
            }
            else
            {
                //if angle_d<0
                //    angle_d=angle_d+180;
                //end
                if (Angle_d < 0)
                {
                    Angle_d += 180;
                }
                //if angle_d>=180
                //    angle_d=angle_d-180;
                //end
                if (Angle_d >= 180)
                {
                    Angle_d -= 180;
                }
            }
            //输出4
            Fi = Angle_d;//fi=angle_d;
            //IPP=[r(2) r(4)];
        }
        static void test_rgpara2()
        {
            //test函数rgpara2
            double[] arrayR = new double[] { -3.1829231e+01, 5.8633108e+00, -3.3810352e+00, -3.8530371e+01 };
            //声明输出
            double Long_a;
            double Short_b;
            double E;
            double Fi;
            rgpara2(arrayR, out Long_a, out Short_b, out E, out Fi);

            Console.WriteLine("计算结果如下: ");
            Console.WriteLine("Long_a = {0}", Long_a);
            Console.WriteLine("Short_b = {0}", Short_b);
            Console.WriteLine("E = {0}", E);
            Console.WriteLine("Fi = {0}", Fi);
            Console.ReadKey();
        }
    }
}
