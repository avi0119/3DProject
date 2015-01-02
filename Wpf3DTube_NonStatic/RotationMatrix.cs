using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Media.Media3D;
namespace Wpf3DTube_NonStatic
{
    interface IMatrixToApplyTRoPoint
    {
        void Calculate3DMatrix(double AngleToRotateBy, double AngleToRotateBy_self, double Tilt, double TiltRelatedRadius);
    }
    public class RotationMatrix : IMatrixToApplyTRoPoint
    {
        // To transform point P by A and then by B, append B to A:
        /*
            C = A.Append(B);
            P' = C.Transform(P);
          */
        protected const double PiConstant = Math.PI;
        private Matrix3D prRotationMatrix3D;
        private Point3D center;
        private Vector3D rotationaxis;
        private double prSinRotorOverCosRotor;
        private double prAngleMovement;
        private double prAngleMovement_self=2;


        public double AngleMovement_self
        {
            get
            {
                return prAngleMovement_self;
            }
            set
            {

                prAngleMovement_self = value;
            }
        }
        public double AngleMovement
        {
            get
            {
                return prAngleMovement;
            }
            set
            {

                prAngleMovement = value;
            }
        }


        public Matrix3D RotationMatrix3D
        {
            get
            {
                return prRotationMatrix3D;
            }
            set
            {

                prRotationMatrix3D = value;
            }
        }

        public Point3D Center
        {
            get
            {
                return center;
            }
            set
            {

                center = value;
            }
        }

        public Vector3D Axis
        {
            get
            {
                return rotationaxis;
            }
            set
            {

                rotationaxis = value / ((value.Length ==0)?1: value.Length);
            }
        }
        public double SinOverCosRotorRatio
        {
            get
            {
                return prSinRotorOverCosRotor;
            }
            set
            {

                prSinRotorOverCosRotor = value;
            }
        }

        //public Matrix3D TheMatrix
        //{
        //    get
        //    {
        //        return prRotationMatrix3D;
        //    }
        //    set
        //    {

        //        prRotationMatrix3D = value;
        //    }
        //}

        public RotationMatrix()
        {

        }
        virtual public void Calculate3DMatrix(double AngleToRotateBy, double AngleToRotateBy_self = 0, double tilt = 0, double prTiltRelatedRadius=0)
        {
            RotationMatrix3D = ReturnRotationMatrix(AngleToRotateBy);
            AngleMovement = AngleToRotateBy;
        }
        public RotationMatrix(Vector3D axis, Point3D CenterPoint=new Point3D(), double SinRotorOverCosRotor = 1)
        {
            Center = CenterPoint;
            Axis = axis;
            SinOverCosRotorRatio = SinRotorOverCosRotor;
            
        }
        virtual public Matrix3D ReturnRotationMatrix(double AngleToRotateBy, bool doitwpfway = false)
        {
            //Matrix3D FinalMatrix;
            //Matrix3D MoveCenterToOrigin;
            //Matrix3D MoveCenterToItsOriginalPosition;
            //Matrix3D RotateAroundAxixComingOutOfOrigin;

            //MoveCenterToOrigin = ReturnMatrix3DToOffsetByXYandZ(-Center.X, -Center.Y, -Center.Z);

            //MoveCenterToItsOriginalPosition = ReturnMatrix3DToOffsetByXYandZ(Center.X, Center.Y, Center.Z);

            //RotateAroundAxixComingOutOfOrigin = ReturnRotationMatrixAroundAxixComingOutOfOrigin(AngleToRotateBy);
            //FinalMatrix = MoveCenterToOrigin;
            //FinalMatrix.Append(RotateAroundAxixComingOutOfOrigin);
            //FinalMatrix.Append(MoveCenterToItsOriginalPosition);
            Vector3D AxisToSpinAround;
            AxisToSpinAround = Axis;
            //AxisToSpinAround = new Vector3D(0, 1, 0);
            Matrix3D ret = ReturnRotationMatrix_inaction(AngleToRotateBy, Center, AxisToSpinAround);
            return ret;
        }
        virtual public Matrix3D ReturnRotationMatrix_inaction(double AngleToRotateBy,Point3D theCenter,Vector3D anAxis,bool clockwise=true)
        {
            Matrix3D FinalMatrix;
            Matrix3D MoveCenterToOrigin;
            Matrix3D MoveCenterToItsOriginalPosition;
            Matrix3D RotateAroundAxixComingOutOfOrigin;

            MoveCenterToOrigin = ReturnMatrix3DToOffsetByXYandZ(-theCenter.X, -theCenter.Y, -theCenter.Z);

            MoveCenterToItsOriginalPosition = ReturnMatrix3DToOffsetByXYandZ(theCenter.X, theCenter.Y, theCenter.Z);

            RotateAroundAxixComingOutOfOrigin = ReturnRotationMatrixAroundAxixComingOutOfOrigin_inaction(anAxis, AngleToRotateBy, clockwise);
            FinalMatrix = MoveCenterToOrigin;
            FinalMatrix.Append(RotateAroundAxixComingOutOfOrigin);
            FinalMatrix.Append(MoveCenterToItsOriginalPosition);

            return FinalMatrix;
        }
        protected Matrix3D ReturnRotationMatrixAroundAxixComingOutOfOrigin(double AngleToRotateBy_indegrees, bool clockWiseDirection=true)
        {
            //double angleInRadians = (AngleToRotateBy_indegrees / 180) * PiConstant;
            //Matrix3D A = new Matrix3D();
            //double COS = Math.Cos(angleInRadians);
            //double SIN = Math.Sin(angleInRadians);
            //double Ux = Axis.X;
            //double Uy = Axis.Y;
            //double Uz= Axis.Z;

            //double Ux2= Math.Pow(Ux,2);
            //double Uy2 = Math.Pow(Uy, 2);
            //double Uz2 = Math.Pow(Uz, 2);
            //double K=SinOverCosRotorRatio;
            //SIN = SIN * K;// in case K is not 1 the rotation is elliptical rather than being circular

            //A.M11 = COS + Ux2*(1-COS);
            //A.M12 = Ux*Uy*(1-COS)-Uz*SIN;
            //A.M13 = Ux * Uz * (1 - COS) + Uy * SIN;
            //A.M14 = 0;

            //A.M21 = Uy*Ux*(1-COS)+Uz*SIN;
            //A.M22 = COS+Uy2*(1-COS);
            //A.M23 = Uy*Uz*(1-COS)-Ux*SIN;
            //A.M24 = 0.0;


            //A.M31 = Uz*Ux*(1-COS)-Uy*SIN;
            //A.M32 = Uz*Uy*(1-COS)+Ux*SIN;
            //A.M33 = COS+Uz2*(1-COS);
            //A.M33 = 0.0;


            //A.M44 = 1.0;
            Matrix3D ret = ReturnRotationMatrixAroundAxixComingOutOfOrigin_inaction(Axis, AngleToRotateBy_indegrees,clockWiseDirection);
            return ret;
        }
        protected Matrix3D ReturnRotationMatrixAroundAxixComingOutOfOrigin_inaction_bkup(Vector3D theAxis,double AngleToRotateBy_indegrees,bool clockWiseDirection=true )
        {
            double angleInRadians = (AngleToRotateBy_indegrees / 180) * PiConstant;
            Matrix3D A = new Matrix3D();
            double COS = Math.Cos(angleInRadians);
            double SIN = Math.Sin(angleInRadians);
            double Ux = theAxis.X;
            double Uy = theAxis.Y;
            double Uz = theAxis.Z;

            double Ux2 = Math.Pow(Ux, 2);
            double Uy2 = Math.Pow(Uy, 2);
            double Uz2 = Math.Pow(Uz, 2);
            double K = SinOverCosRotorRatio;
            //K = 1;
            SIN = SIN * K;// in case K is not 1 the rotation is elliptical rather than being circular
            if (clockWiseDirection == false) SIN = SIN * -1;
            A.M11 = COS + Ux2 * (1 - COS);
            A.M12 = Ux * Uy * (1 - COS) - Uz * SIN;
            A.M13 = Ux * Uz * (1 - COS) + Uy * SIN;
            A.M14 = 0;

            A.M21 = Uy * Ux * (1 - COS) + Uz * SIN;
            A.M22 = COS + Uy2 * (1 - COS);
            A.M23 = Uy * Uz * (1 - COS) - Ux * SIN;
            A.M24 = 0.0;


            A.M31 = Uz * Ux * (1 - COS) - Uy * SIN;
            A.M32 = Uz * Uy * (1 - COS) + Ux * SIN;
            A.M33 = COS + Uz2 * (1 - COS);
            A.M34 = 0.0;


            A.M44 = 1.0;
            A.OffsetX = A.OffsetY = A.OffsetZ = 0;
            //A = RetMatrixPerRotationAroundAxisGivemCenterAndAngle(-AngleToRotateBy_indegrees, new Point3D(0, 0, 0), theAxis);
            return A;
        }
        protected Matrix3D ReturnRotationMatrixAroundAxixComingOutOfOrigin_inaction(Vector3D theAxis, double AngleToRotateBy_indegrees, bool clockWiseDirection = true)
        {
            double angleInRadians = (AngleToRotateBy_indegrees / 180) * PiConstant;
            Matrix3D A = new Matrix3D();
            double COS = Math.Cos(angleInRadians);
            double SIN = Math.Sin(angleInRadians);
            double Ux = theAxis.X;
            double Uy = theAxis.Y;
            double Uz = theAxis.Z;

            double Ux2 = Math.Pow(Ux, 2);
            double Uy2 = Math.Pow(Uy, 2);
            double Uz2 = Math.Pow(Uz, 2);
            double K = SinOverCosRotorRatio;
            //K = 1;
            SIN = SIN * K;// in case K is not 1 the rotation is elliptical rather than being circular
            if (clockWiseDirection == false) SIN = SIN * -1;
            A.M11 = COS + Ux2 * (1 - COS);
            A.M12 = Ux * Uy * (1 - COS) - Uz * SIN;
            A.M13 = -(Ux * Uz * (1 - COS) + Uy * SIN);//
            A.M14 = 0;

            A.M21 = Uy * Ux * (1 - COS) + Uz * SIN;
            A.M22 = COS + Uy2 * (1 - COS);
            A.M23 = -(Uy * Uz * (1 - COS) - Ux * SIN);//
            A.M24 = 0.0;


            A.M31 = -(Uz * Ux * (1 - COS) - Uy * SIN);//
            A.M32 = -(Uz * Uy * (1 - COS) + Ux * SIN);//
            A.M33 = COS + Uz2 * (1 - COS);
            A.M34 = 0.0;


            A.M44 = 1.0;
            A.OffsetX = A.OffsetY = A.OffsetZ = 0;
            //A = RetMatrixPerRotationAroundAxisGivemCenterAndAngle(-AngleToRotateBy_indegrees, new Point3D(0, 0, 0), theAxis);
            return A;
        }

        public static  Matrix3D ReturnRotationMatrixAroundAxixComingOutOfOrigin_inactionstatic_bkup(Vector3D theAxis, double AngleToRotateBy_indegrees, bool clockWiseDirection = true)
        {
            double angleInRadians = (AngleToRotateBy_indegrees / 180) * PiConstant;
            Matrix3D A = new Matrix3D();
            double COS = Math.Cos(angleInRadians);
            double SIN = Math.Sin(angleInRadians);
            double Ux = theAxis.X;
            double Uy = theAxis.Y;
            double Uz = theAxis.Z;

            double Ux2 = Math.Pow(Ux, 2);
            double Uy2 = Math.Pow(Uy, 2);
            double Uz2 = Math.Pow(Uz, 2);
            double K = 1;
            //K = 1;
            SIN = SIN * K;// in case K is not 1 the rotation is elliptical rather than being circular
            if (clockWiseDirection == false) SIN = SIN * -1;
            A.M11 = COS + Ux2 * (1 - COS);
            A.M12 = Ux * Uy * (1 - COS) - Uz * SIN;
            A.M13 = Ux * Uz * (1 - COS) + Uy * SIN;
            A.M14 = 0;

            A.M21 = Uy * Ux * (1 - COS) + Uz * SIN;
            A.M22 = COS + Uy2 * (1 - COS);
            A.M23 = Uy * Uz * (1 - COS) - Ux * SIN;
            A.M24 = 0.0;


            A.M31 = Uz * Ux * (1 - COS) - Uy * SIN;
            A.M32 = Uz * Uy * (1 - COS) + Ux * SIN;
            A.M33 = COS + Uz2 * (1 - COS);
            A.M34 = 0.0;


            A.M44 = 1.0;
            A.OffsetX = A.OffsetY = A.OffsetZ = 0;
            //A = RetMatrixPerRotationAroundAxisGivemCenterAndAngle(-AngleToRotateBy_indegrees, new Point3D(0, 0, 0), theAxis);
            return A;
        }
        public static Matrix3D ReturnRotationMatrixAroundAxixComingOutOfOrigin_inactionstatic(Vector3D theAxis, double AngleToRotateBy_indegrees, bool clockWiseDirection = true)
        {
            double angleInRadians = (AngleToRotateBy_indegrees / 180) * PiConstant;
            Matrix3D A = new Matrix3D();
            double COS = Math.Cos(angleInRadians);
            double SIN = Math.Sin(angleInRadians);
            double Ux = theAxis.X;
            double Uy = theAxis.Y;
            double Uz = theAxis.Z;

            double Ux2 = Math.Pow(Ux, 2);
            double Uy2 = Math.Pow(Uy, 2);
            double Uz2 = Math.Pow(Uz, 2);
            double K = 1;
            //K = 1;
            SIN = SIN * K;// in case K is not 1 the rotation is elliptical rather than being circular
            if (clockWiseDirection == false) SIN = SIN * -1;
            A.M11 = COS + Ux2 * (1 - COS);
            A.M12 = Ux * Uy * (1 - COS) - Uz * SIN;
            A.M13 =-( Ux * Uz * (1 - COS) + Uy * SIN);//
            A.M14 = 0;

            A.M21 = Uy * Ux * (1 - COS) + Uz * SIN;
            A.M22 = COS + Uy2 * (1 - COS);
            A.M23 =-( Uy * Uz * (1 - COS) - Ux * SIN);//
            A.M24 = 0.0;


            A.M31 = -(Uz * Ux * (1 - COS) - Uy * SIN);//
            A.M32 = -(Uz * Uy * (1 - COS) + Ux * SIN);//
            A.M33 = COS + Uz2 * (1 - COS);
            A.M34 = 0.0;


            A.M44 = 1.0;
            A.OffsetX = A.OffsetY = A.OffsetZ = 0;
            //A = RetMatrixPerRotationAroundAxisGivemCenterAndAngle(-AngleToRotateBy_indegrees, new Point3D(0, 0, 0), theAxis);
            return A;
        }
        private Matrix3D RetMatrixPerRotationAroundAxisGivemCenterAndAngle(double AngleToRotateBy, Point3D aCenter, Vector3D anAxis)
        {
            Rotation3D b = new AxisAngleRotation3D(anAxis, AngleToRotateBy);
            RotateTransform3D a = new RotateTransform3D(b, aCenter);
            return a.Value;

        }

        protected Matrix3D ReturnMatrix3DToOffsetByXYandZ(double x, double y, double z)
        {
            Matrix3D m = new Matrix3D();
            m.M11 = 1.0; m.M12 = 0.0; m.M13 = 0.0; m.M14 = 0.0;
            m.M21 = 0.0; m.M22 = 1.0; m.M23 = 0.0; m.M24 = 0.0;
            m.M31 = 0.0; m.M32 = 0.0; m.M33 = 1.0; m.M34 = 0.0;
            m.OffsetX = x; m.OffsetY = y; m.OffsetZ = z; m.M44 = 1.0;

            return m;
        }
        public static Matrix3D ReturnMatrix3DToOffsetByXYandZ_static(double x, double y, double z)
        {
            Matrix3D m = new Matrix3D();
            m.M11 = 1.0; m.M12 = 0.0; m.M13 = 0.0; m.M14 = 0.0;
            m.M21 = 0.0; m.M22 = 1.0; m.M23 = 0.0; m.M24 = 0.0;
            m.M31 = 0.0; m.M32 = 0.0; m.M33 = 1.0; m.M34 = 0.0;
            m.OffsetX = x; m.OffsetY = y; m.OffsetZ = z; m.M44 = 1.0;

            return m;
        }
        protected Matrix3D ReturnMatrix3DToOffsetByXYandZ_old(double x, double y, double z)
        {
            Matrix3D A = new Matrix3D();
            A.M14 = 0;
            A.M24 = 0;
            A.M34 = 0;
            A.OffsetX = x;
            A.OffsetY = y;
            A.OffsetZ = z;


            A.M11 = 1;
            A.M22 = 1;
            A.M33 = 1;
            A.M44 = 1.0;

            return A;
        }
        

    }
    public class RotationMatrix_v2 : RotationMatrix
    {
        protected Vector3D prTilitAxix;
        public Vector3D TilitAxix
        {
            get
            {
                return prTilitAxix;
            }
            set
            {

                prTilitAxix = value;
            }
        }

        protected bool prTilitAxixFlag;
        public bool TilitAxixFlag
        {
            get
            {
                return prTilitAxixFlag;
            }
            set
            {

                prTilitAxixFlag = value;
            }
        }
        public Dictionary<string, Point3DOriginalAndNow> Important3DPoints;

        public RotationMatrix_v2()
        {
        }
        public RotationMatrix_v2(Vector3D axis, Point3D CenterPoint = new Point3D(), double SinRotorOverCosRotor = 1):base( axis,  CenterPoint ,  SinRotorOverCosRotor)
        {
        }
        override public void Calculate3DMatrix(double AngleToRotateBy, double AngleToRotateBy_self, double tilt, double TiltRelatedRadius)
        {

            Point3D p = new Point3D(0, 1, 0);
            Point3D result;
            Matrix3D RotationMatrixToApplyOnCenterPointOnly = ReturnRotationMatrix(AngleToRotateBy);
            //Point3D SatteliteCenterPoint=f
            foreach (string PointName in Important3DPoints.Keys)
            {
                Point3DOriginalAndNow CurrentOriginalAndNow = Important3DPoints[PointName];
                CurrentOriginalAndNow.Now = RotationMatrixToApplyOnCenterPointOnly.Transform(CurrentOriginalAndNow.Original);

            }
            Vector3D offsetVector = (Vector3D)Important3DPoints["CENTER POINT"].Offset;
            Random random = new Random();
            //double xmoverandom = random.NextDouble()*100;
            //double ymoverandom = random.NextDouble() * 100;
            //double zmoverandom = random.NextDouble() * 100;
            ////offsetVector.X = xmoverandom;
            ////offsetVector.Y = ymoverandom;
            ////offsetVector.Z = zmoverandom;
            //Matrix3D FinalMatrix;
            Matrix3D MoveCenterToOrigin;


            MoveCenterToOrigin = ReturnMatrix3DToOffsetByXYandZ(offsetVector.X, offsetVector.Y, offsetVector.Z);
            result = MoveCenterToOrigin.Transform(p);
            result = result;
            Matrix3D FinalMatrix;
            Matrix3D TiltMAtrix;
            Matrix3D RotateAroundSelfAxis;
            FinalMatrix = MoveCenterToOrigin;
            //AngleToRotateBy_self
            RotateAroundSelfAxis = ReturnRotationMatrix_inaction(AngleToRotateBy_self, Important3DPoints["CENTER POINT"].Now, Axis, true);
            //Point3D a = new Point3D(10, 0, 0);
            //Point3D b= RotateAroundSelfAxis.Transform(a);
            //Point3D c = b;
            //if (true == false)
            //{
            //this matrix is responsible for self rotation
            FinalMatrix.Append(RotateAroundSelfAxis);
            //}
            //FinalMatrix.Append(RotateAroundSelfAxis);
            double radius = TiltRelatedRadius;
            TiltMAtrix = TiltGivenAxisByXDegrees(Axis, Important3DPoints["CENTER POINT"].Now, tilt, radius);
            TiltMAtrix = ReturnAxisToTiltShapeAround(Axis, Important3DPoints["CENTER POINT"].Now, tilt, radius);
            if (!(tilt == 0))
            {
                FinalMatrix.Append(TiltMAtrix);
            }
            RotationMatrix3D = FinalMatrix;
            RotationMatrix3D = ApplyAllMatrices(MoveCenterToOrigin, RotateAroundSelfAxis, TiltMAtrix, tilt); ;

        }
        private Matrix3D ApplyAllMatrices(Matrix3D MatrixTranslation, Matrix3D MatrixSelfRotation, Matrix3D MatrixTilt, double tilt) 
        {
            Matrix3D identiry = Matrix3D.Identity;
            Matrix3D[] matrixArray = new Matrix3D [3];//{ MatrixTranslation, MatrixSelfRotation, tilt == 0 ? identiry : MatrixTilt };
            //Matrix3D[] matrixArray = new Matrix3D { MatrixTranslation, MatrixSelfRotation, tilt==0 ? identiry:MatrixTilt };
            matrixArray[0] =  MatrixTranslation ;
            matrixArray[1] = identiry;// MatrixSelfRotation;
            matrixArray[2] =   tilt == 0 ? identiry : MatrixTilt ;
            string order = "1,2,3";
            Matrix3D FinalMatrix;
            string[] words = order.Split(',');
            int[] ia = words.Select(n => Convert.ToInt32(n)).ToArray();
            FinalMatrix = matrixArray[(ia[0]-1)];
            for (int counter = 1; counter < ia.Length; ++counter)
            {
                FinalMatrix.Append(matrixArray[(ia[counter]-1)]);
            }
            return FinalMatrix;
        }

         public void Calculate3DMatrix_bkup2(double AngleToRotateBy, double AngleToRotateBy_self, double tilt, double TiltRelatedRadius)
        {

            Point3D p = new Point3D(0, 1, 0);
            Point3D result;
            Matrix3D RotationMatrixToApplyOnCenterPointOnly = ReturnRotationMatrix(AngleToRotateBy);
            //Point3D SatteliteCenterPoint=f
            foreach (string PointName in Important3DPoints.Keys)
            {
                Point3DOriginalAndNow CurrentOriginalAndNow = Important3DPoints[PointName];
                CurrentOriginalAndNow.Now = RotationMatrixToApplyOnCenterPointOnly.Transform(CurrentOriginalAndNow.Original);

            }
            Vector3D offsetVector = (Vector3D)Important3DPoints["CENTER POINT"].Offset;

            //Matrix3D FinalMatrix;
            Matrix3D MoveCenterToOrigin;


            MoveCenterToOrigin = ReturnMatrix3DToOffsetByXYandZ(offsetVector.X, offsetVector.Y, offsetVector.Z);
            result = MoveCenterToOrigin.Transform(p);
            result = result;
            Matrix3D FinalMatrix;
            Matrix3D TiltMAtrix;
            Matrix3D RotateAroundSelfAxis;
            FinalMatrix = MoveCenterToOrigin;
            //AngleToRotateBy_self
            RotateAroundSelfAxis = ReturnRotationMatrix_inaction(AngleToRotateBy_self, Important3DPoints["CENTER POINT"].Now, Axis, true);
            //Point3D a = new Point3D(10, 0, 0);
            //Point3D b= RotateAroundSelfAxis.Transform(a);
            //Point3D c = b;
            //if (true == false)
            //{
            //this matrix is responsible for self rotation
            FinalMatrix.Append(RotateAroundSelfAxis);
            //}
            //FinalMatrix.Append(RotateAroundSelfAxis);
            double radius = TiltRelatedRadius;
            TiltMAtrix = TiltGivenAxisByXDegrees(Axis, Important3DPoints["CENTER POINT"].Now, tilt, radius);
            TiltMAtrix = ReturnAxisToTiltShapeAround(Axis, Important3DPoints["CENTER POINT"].Now, tilt, radius);
            if (!(tilt == 0))
            {
                FinalMatrix.Append(TiltMAtrix);
            }
            RotationMatrix3D = FinalMatrix;

        }
        public void Calculate3DMatrix_bkup(double AngleToRotateBy, double AngleToRotateBy_self, double tilt, double TiltRelatedRadius)

        {

            Point3D p = new Point3D(0, 1, 0);
            Point3D result;
            Matrix3D RotationMatrixToApplyOnCenterPointOnly = ReturnRotationMatrix(AngleToRotateBy);
            //Point3D SatteliteCenterPoint=f
            foreach (string PointName in Important3DPoints.Keys)
            {
                Point3DOriginalAndNow CurrentOriginalAndNow = Important3DPoints[PointName];
                CurrentOriginalAndNow.Now = RotationMatrixToApplyOnCenterPointOnly.Transform(CurrentOriginalAndNow.Original);

            }
            Vector3D offsetVector = (Vector3D)Important3DPoints["CENTER POINT"].Offset;

            //Matrix3D FinalMatrix;
            Matrix3D MoveCenterToOrigin;


            MoveCenterToOrigin = ReturnMatrix3DToOffsetByXYandZ(offsetVector.X, offsetVector.Y, offsetVector.Z);
            result = MoveCenterToOrigin.Transform(p);
            result = result;
            Matrix3D FinalMatrix;
            Matrix3D TiltMAtrix;
            Matrix3D RotateAroundSelfAxis;
            FinalMatrix = MoveCenterToOrigin;
                                                                  //AngleToRotateBy_self
            RotateAroundSelfAxis = ReturnRotationMatrix_inaction(AngleToRotateBy_self, Important3DPoints["CENTER POINT"].Now, Axis, true);
            //Point3D a = new Point3D(10, 0, 0);
            //Point3D b= RotateAroundSelfAxis.Transform(a);
            //Point3D c = b;
            FinalMatrix.Append(RotateAroundSelfAxis);
            //FinalMatrix.Append(RotateAroundSelfAxis);
            double radius = TiltRelatedRadius;
            TiltMAtrix = TiltGivenAxisByXDegrees(Axis, Important3DPoints["CENTER POINT"].Now, tilt, radius);
            TiltMAtrix = ReturnAxisToTiltShapeAround(Axis, Important3DPoints["CENTER POINT"].Now, tilt, radius);
            if (!(tilt==0)) {
                FinalMatrix.Append(TiltMAtrix);
            }
            RotationMatrix3D = FinalMatrix;

        }
        private Matrix3D TiltGivenAxisByXDegrees(Vector3D anAxisPerpendicular, Point3D currentCenterPoint, double NumberOfDegreesToTiltBy,double radius) 
        {

            /// radius is half of the perpendicular axis.  we use it to go back to the beginning of the axis
            /// and then tilt it from that same beginning. When it comes to eliipsoid, the beginning is the center minus
            /// // "half of the perpendicular"
            anAxisPerpendicular = anAxisPerpendicular / anAxisPerpendicular.Length;
            Point3D currentCenterPoint2 = currentCenterPoint - anAxisPerpendicular * radius;
            Vector3D RotorFromBAseOfAxisToOrigin = (Vector3D)(currentCenterPoint2) * -1;
            double angleInRadians = (NumberOfDegreesToTiltBy / 180) * PiConstant;

            double COS = Math.Cos(angleInRadians);
            double SIN = Math.Sin(angleInRadians);
            
            RotorFromBAseOfAxisToOrigin = RotorFromBAseOfAxisToOrigin / RotorFromBAseOfAxisToOrigin.Length;
            Vector3D NewAxis = anAxisPerpendicular * COS + RotorFromBAseOfAxisToOrigin*SIN;
            //return ret/ret.Length; 
            Point3D newcenter = currentCenterPoint2 + NewAxis / NewAxis.Length * radius;
            Vector3D offset = (Vector3D)newcenter - (Vector3D)Important3DPoints["CENTER POINT"].Now;
            Matrix3D FinalMatrix;
            FinalMatrix = ReturnMatrix3DToOffsetByXYandZ(offset.X, offset.Y, offset.Z);
            return FinalMatrix;
        }
        private Matrix3D ReturnAxisToTiltShapeAround( Vector3D anAxisPerpendicular, Point3D currentCenterPoint, double NumberOfDegreesToTiltBy, double radius)
        {
            
  

            anAxisPerpendicular = anAxisPerpendicular / anAxisPerpendicular.Length;
            Point3D currentCenterPoint2 = currentCenterPoint - anAxisPerpendicular * radius;
            Point3D CenterOfReturnAxis = (Point3D)currentCenterPoint2;
            if (TilitAxixFlag==false)
            {
                Vector3D retAxis;
                Vector3D RotorFromBAseOfAxisToOrigin = (Vector3D)(currentCenterPoint);
                RotorFromBAseOfAxisToOrigin = RotorFromBAseOfAxisToOrigin / RotorFromBAseOfAxisToOrigin.Length;
                retAxis = Vector3D.CrossProduct(anAxisPerpendicular, RotorFromBAseOfAxisToOrigin);

                TilitAxix = retAxis / retAxis.Length;
                TilitAxixFlag = true;
            }
            Matrix3D FinalMatrix = ReturnRotationMatrix_inaction(NumberOfDegreesToTiltBy, CenterOfReturnAxis, TilitAxix); ;
            //return retAxis / retAxis.Length;
            //double angleInRadians = (NumberOfDegreesToTiltBy / 180) * PiConstant;
            
            //double COS = Math.Cos(angleInRadians);
            //double SIN = Math.Sin(angleInRadians);

            
            //Vector3D NewAxis = anAxisPerpendicular * COS + RotorFromBAseOfAxisToOrigin * SIN;
            //return ret/ret.Length; 
            //Point3D newcenter = currentCenterPoint2 + NewAxis / NewAxis.Length * radius;
            //Vector3D offset = (Vector3D)newcenter - (Vector3D)Important3DPoints["CENTER POINT"].Now;
            //Matrix3D FinalMatrix;
            //FinalMatrix = ReturnMatrix3DToOffsetByXYandZ(offset.X, offset.Y, offset.Z);
            return FinalMatrix;
        }

        override public Matrix3D ReturnRotationMatrix(double AngleToRotateBy,bool doitwpfway=false)
        {
            Matrix3D FinalMatrix_old;
            Matrix3D FinalMatrix;
            Matrix3D MoveCenterToOrigin;
            Matrix3D MoveCenterToItsOriginalPosition;
            Matrix3D RotateAroundAxixComingOutOfOrigin;
  
                MoveCenterToOrigin = ReturnMatrix3DToOffsetByXYandZ(-Center.X, -Center.Y, -Center.Z);

                MoveCenterToItsOriginalPosition = ReturnMatrix3DToOffsetByXYandZ(Center.X, Center.Y, Center.Z);

                RotateAroundAxixComingOutOfOrigin = ReturnRotationMatrixAroundAxixComingOutOfOrigin(AngleToRotateBy);
                FinalMatrix_old = MoveCenterToOrigin;
                FinalMatrix_old.Append(RotateAroundAxixComingOutOfOrigin);
                FinalMatrix_old.Append(MoveCenterToItsOriginalPosition);

                if (doitwpfway == true)
                {
                    Rotation3D b = new AxisAngleRotation3D(Axis, AngleToRotateBy);
                    RotateTransform3D a = new RotateTransform3D(b, Center);

                    FinalMatrix = a.Value;
                }
                else
                {
                    FinalMatrix = FinalMatrix_old;
                }

            return FinalMatrix;
        }


    }
}
