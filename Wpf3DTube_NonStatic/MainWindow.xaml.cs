#define   USEVERSION2
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
using System.Drawing;
using System.Resources;


namespace Wpf3DTube_NonStatic
{
    
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    
    public static class GlobalVar
    {
        public const bool USEVERSION2 = true;
        public const    double GravitationalG=6.67384 * 10-11;// m3 kg-1 s-2;
        public static MainWindow mainwindow;
        static public Vector3D CalculateTangerntenialVector(Vector3D PrePos, Vector3D CurrentPos,double  NumberOfUnitOfTiems)
        {
            Vector3D a = CurrentPos - PrePos;
            return a / NumberOfUnitOfTiems;
        }
    }
    public class OrbittingPlane
    {
        public PlanetPropeties planet { get; set; }// based on normal of planet
        public double inclinationWithRespectToNormal { get; set; }

    }
    public class PlanetDimensions
    {
        public double EquaterialSemiMajorLength { get; set; }
        public double EquaterialSemiMinorLength { get; set; }
        public double PolarRadius { get; set; }


    }
    public class OrbitAroundData
    {
        public PlanetPropeties planet { get; set; }
        public Plane intersectionPlane { get; set; }
        public bool usePlanetLocatonAsPointInIntersectingLine { get; set; }
        public OrbittingPlane OrbittingPlaneAngleDeterminator { get; set; }
    }
    public class PlanetPropeties
    {
        public const double VisualToRealRAtio=10.0;
        public PlanetDimensions Dimensions  { get; set; } 
        public double aphelion { get; set; } //max distance
        public double perihelion { get; set; } //min distance
        public double mass { get; set; } //min distance
        public double SpeedAtperihelion { get; set; } //MAX SPEED
        public double SpeedAtaphelion { get; set; } //Min SPEED
        public double OrbitSemiMinorAxis { get; set; }
        public double OrbitSemiMajorAxis { get; set; }
        public Vector3D SemiMajorAxixProjectionOntoReferncePlane { get; set; }
        public Vector3D SemiMajorAxisDirection { get; set; }
        public Vector3D SemiMinorAxisDirection { get; set; }
        public double SinRotorOverCosRotor { get; set; }
        public double Tilt { get; set; }
        public double eccentricity { get; set; }
        public Point3D imaginaryCenter { get; set; }
        public Point3D otherFocalPoint { get; set; }
        public Point3D initilaLocation { get; set; }
        public Point3D myCurrentLocation { get; set; }
        public double angleWithRespectToInitialLocation { get; set; }// refers to initioal location which is assumed to be on semimajor axix
        public Vector3D mormalToPlanetEquator { get; set; }

        public OrbitAroundData MyOrbitFocalPlanet { get; set; }
        public Plane PlaneIntersectingOrbit { get; set; }
        public Point3D calcultePositionAtGivenAngle(double angleWithRespectToSemiMajor)
        {

            double angleInRadians = angleWithRespectToSemiMajor / 180 * Math.PI;
            double Cos = Math.Cos(angleInRadians);
            double Sin = Math.Sin(angleInRadians);

            Vector3D C = SemiMajorAxisDirection; Vector3D S = SemiMinorAxisDirection;
            C.Normalize(); 
            S.Normalize();
            Point3D p = imaginaryCenter + Cos * C + Sin * S;

            return p;
        }
        public Vector3D returnTangetToGivenPointOnCircumference(Point3D p) 
        {
            Point3D f1 = MyOrbitFocalPlanet.planet.myCurrentLocation;
            Point3D f2 = otherFocalPoint;
            Vector3D dist1 = f1-p;
            Vector3D dist2 = f2-p;
            Vector3D lineconnectingf1andf2 = f1 - f2;
            double lenght_dist1 = dist1.Length;
            double lenght_dist2 = dist2.Length;
            double length_lineconnectingf1andf2 = lineconnectingf1andf2.Length;
            double anglebtewwndist1anddist2_part1 = Math.Pow(length_lineconnectingf1andf2, 2) - Math.Pow(lenght_dist1, 2) - Math.Pow(lenght_dist2, 2);
            double anglebtewwndist1anddist2_part12 = -2 * lenght_dist1 * lenght_dist2;
            double anglebtewwndist1anddist2_part13 = anglebtewwndist1anddist2_part1 / anglebtewwndist1anddist2_part12;
            double nglebtewwndist1anddist2_cosInRadians=Math.Acos(anglebtewwndist1anddist2_part13);
            double angleInDegrees = nglebtewwndist1anddist2_cosInRadians / Math.PI * 180;
            double angleBetweenTangentAmdDist1or2 = (180 - angleInDegrees) / 2;
            GramSchmidt GM = new GramSchmidt(true, dist1, dist2);
            Vector3D[] basis_orthonormalOFDists=(GM.OrthonormalBasis).ToArray();// one of them is one of the dists, the other is perpendicular. we will need to rotate them by angleBetweenTangentAmdDist1or2

            double Cos = Math.Cos(angleBetweenTangentAmdDist1or2);
            double Sin = Math.Sin(angleBetweenTangentAmdDist1or2);
            Vector3D v = basis_orthonormalOFDists[0]*Cos + basis_orthonormalOFDists[1]*Sin;
            return v;
        }
        public void calculateSemiMinor()
        {
            OrbitSemiMinorAxis = OrbitSemiMajorAxis * Math.Sqrt(1 - Math.Pow(eccentricity, 2));
            SinRotorOverCosRotor = OrbitSemiMinorAxis / OrbitSemiMajorAxis;
        }
        public void CalcualteOtherFocalPoint()
        {
            double focalPointDistanceFromCenter = Math.Sqrt((Math.Pow(OrbitSemiMajorAxis, 2) - Math.Pow(OrbitSemiMinorAxis, 2)));
            otherFocalPoint = (Point3D)(imaginaryCenter + SemiMajorAxisDirection * focalPointDistanceFromCenter);
        }
        public void CalcualteImaginaryCenter()
        {
            imaginaryCenter = (Point3D)(MyOrbitFocalPlanet.planet.initilaLocation - perihelion * SemiMajorAxisDirection + SemiMajorAxisDirection * OrbitSemiMajorAxis);
            
        }
        public void CalculateSemiMajorAxixProjectionOntoReferncePlane()
        {
            Plane aPlane = MyOrbitFocalPlanet.intersectionPlane;
            SemiMajorAxixProjectionOntoReferncePlane = (Plane.ReturnIntersectionLine(aPlane, new Plane(MyOrbitFocalPlanet.planet.mormalToPlanetEquator, MyOrbitFocalPlanet.planet.initilaLocation))).VectorDirection;
            SemiMajorAxixProjectionOntoReferncePlane = SemiMajorAxixProjectionOntoReferncePlane / SemiMajorAxixProjectionOntoReferncePlane.Length;
        
        }
        public void calculateNormalToEquatorAndSemiMAjorAxisDirection( )
        {
            /// this is before tilt is applied
            /// 
            PlanetPropeties ReferncePlanet = this.MyOrbitFocalPlanet.OrbittingPlaneAngleDeterminator.planet;
            double inclination = this.MyOrbitFocalPlanet.OrbittingPlaneAngleDeterminator.inclinationWithRespectToNormal;
            Plane planeOfReferncePlanet = new Plane(ReferncePlanet.mormalToPlanetEquator, ReferncePlanet.initilaLocation);
            Vector3D MyMajorAxisVecorProjection = SemiMajorAxixProjectionOntoReferncePlane;// returnVecotorWithWhichMajorAxisCoincides();// this.OrbitSemiMajorAxis;
            Vector3D slantedVector = ReturnVectorWithXAngleBetweenSaidVectorAndPlane(inclination, MyMajorAxisVecorProjection, ReferncePlanet.mormalToPlanetEquator);
            SemiMajorAxisDirection = slantedVector / slantedVector.Length;
            Vector3D OrtihgonalToProjecteDaNDsLANTED = Vector3D.CrossProduct(slantedVector, MyMajorAxisVecorProjection);
            //OrtihgonalToProjecteDaNDsLANTED=OrtihgonalToProjecteDaNDsLANTED/OrtihgonalToProjecteDaNDsLANTED.Length;
            Vector3D requiredNormal;
            if (OrtihgonalToProjecteDaNDsLANTED.Length == 0)
            {
                requiredNormal=-ReferncePlanet.mormalToPlanetEquator;
            }
            else
            {
                //our sought vector is that one which is perpendiculat to the slanted and that which is orthogonal to the slanted and the semimahor projection
                requiredNormal = Vector3D.CrossProduct(OrtihgonalToProjecteDaNDsLANTED, slantedVector);
            }
            this.mormalToPlanetEquator =requiredNormal/-requiredNormal.Length;
            SemiMinorAxisDirection = Vector3D.CrossProduct(this.mormalToPlanetEquator, SemiMajorAxisDirection);
            SemiMinorAxisDirection.Normalize();

        }
        public Vector3D returnVecotorWithWhichMajorAxisCoincides()
        {
            /// vector between the centers of earth, for example, and the center of sun which is the focal point
            Vector3D v = imaginaryCenter - this.MyOrbitFocalPlanet.planet.imaginaryCenter;

            return v;
        }
        private Vector3D ReturnVectorWithXAngleBetweenSaidVectorAndPlane(double x, Vector3D vectorToWhichInclinationRefersTo,Vector3D normal)
        {
            //Vector3D vectorPerPendicularToBothGivenVectors = Vector3D.CrossProduct(vectorToWhichInclinationRefersTo, normal);
            //vectorPerPendicularToBothGivenVectors = vectorPerPendicularToBothGivenVectors/vectorPerPendicularToBothGivenVectors.Length;
            //// we have the sun's  normal to its equator plane plus a vector which is on that plane (preoendicular to the normal) and which is a projection
            /// of the major axis vecotr with magnitude of 1
            // now we find the vector which will be on the plane of the body in question taking into axccount the inclination
            double PI = Math.PI;
            double sinA = Math.Sin(x/180*PI);
            double cosA = Math.Cos(x / 180 * PI);
            // weassume the normal and the projection vector are unti vectors;
            Vector3D v1 = normal * sinA;
            Vector3D v2 = vectorToWhichInclinationRefersTo * cosA;
            return (v1 + v2)/(v1+v2).Length;
            
        }

    }
    public class SHapeParameters
    {
        //private const PlanetPropeties constPlanetPropeties = new PlanetPropeties();
        public bool prDisplayCommet;
        public Vector3D prAxisVector;
        public bool prIsItTile;
        public Vector3D prTranslationVector;
        public double prTubeThickness;
        public double prInternalDiameter;
        public int prangleStep_horizotal;
        public int prangleStep_vertical;
        public int prNumberOfDefreesToComplete_horizontal;
        public int prNumberOfDefreesToComplete_vertical;
        public object prImageBrush;
        public double prinitialselfRotation;
        public double prinitialGroundAngle;
        public Type prShapeType;
        public double prEllipsisA;
        public double prEllipsisB;
        public double prEllipsisC;
        public double prAngleMovementAroundExternalAxis;
        public double prAngleMovementAroundInternalAxis;
        public double prDistanceFromRespectiveSun;
        public double prInitialOrbitalAngle;
        public double prTilt;
        public bool prDisplayAxixPole;
        public bool prDisplayTrail;
        public PlanetPropeties prPlanetRelated;


        public SHapeParameters(
            Vector3D AxisVector,
            Vector3D TranslationVector,
            double TubeThickness,
            double InternalDiameter,
            int angleStep_horizotal,
            int angleStep_vertical,
            int NumberOfDefreesToComplete_horizontal,
            int NumberOfDefreesToComplete_vertical,
            object aImageBrush,
            double initialselfRotation,
            double initialGroundAngle,
            bool rIsItTile,
            bool DisplayCommet,
            Type ShapeType,
            PlanetPropeties PlanetRelated,
            double EllipsisA=0,
            double EllipsisB=0,
            double EllipsisC=0,
            double AngleMovementAroundExternalAxis=1,
            double AngleMovementAroundInternalAxis=1,
            double DistanceFromRespectiveSun=10,
            double InitialOrbitalAngle = 0,
            double Tilt=0,
            bool DisplayAxixPole=false,
            bool DisplayTrail = false
            
            )
        {
            prAxisVector = AxisVector;
            prTranslationVector = TranslationVector;
            prTubeThickness = TubeThickness;
            prInternalDiameter = InternalDiameter;
            prangleStep_horizotal = angleStep_horizotal;
            prangleStep_vertical = angleStep_vertical;
            prNumberOfDefreesToComplete_horizontal = NumberOfDefreesToComplete_horizontal;
            prNumberOfDefreesToComplete_vertical = NumberOfDefreesToComplete_vertical;
            prImageBrush = aImageBrush;
            prinitialselfRotation=initialselfRotation;
            prinitialGroundAngle = initialGroundAngle;
            prIsItTile = rIsItTile;
            prDisplayCommet = DisplayCommet;
            prShapeType = ShapeType;
            prEllipsisA = EllipsisA;
            prEllipsisB = EllipsisB;
            prEllipsisC = EllipsisC;
            prAngleMovementAroundExternalAxis = AngleMovementAroundExternalAxis;
            prAngleMovementAroundInternalAxis = AngleMovementAroundInternalAxis;
            prDistanceFromRespectiveSun = DistanceFromRespectiveSun;
            prInitialOrbitalAngle = InitialOrbitalAngle;
            prTilt = Tilt;
            prDisplayAxixPole = DisplayAxixPole;
            prDisplayTrail = DisplayTrail;
            prPlanetRelated = PlanetRelated;


        }
    }
    public partial class MainWindow : Window
    {
        public static double angleaccumalated_delet=0;
        //public AxisAngleRotation3D rotationY = new AxisAngleRotation3D();
        //public MatrixTransform3D myMatrixTransform3D;
        //static private double selfrotatingangle = 0;
        //static private double groundrotatingangle = 0;
        //TubeGeometry TG;
        public const double AstronomicalUnint = 150000;
        Dictionary<string, PlanetPropeties> Planets = new Dictionary<string, PlanetPropeties>();
        static List<TubeGeometry> ShapesList = new List<TubeGeometry>();
        private GeometryModel3D mGeometry;
        private GeometryModel3D mGeometry2; 
        private GeometryModel3D mGeometry3;
        static TimeSpan lastRenderTime = new TimeSpan();
        //public static OrthographicCamera Mycamera;
        public static PerspectiveCamera Mycamera;
        private int SpinMode = 0;
        private Model3DGroup arrowbodyandhead;
        private Model3DGroup CubeFrame;
        static public  Model3DGroup globalTrialCube;
        static public MeshGeometry3D globalMesh;
        static public  Viewport3D globalviewport;
        static public bool displayRay;
        public Visual3D lookDirectionRay;
        public MainWindow()
        {
            displayRay = false;
            testPlanes();
            //testGramSchmidt();
            angleaccumalated_delet = 0;
            GlobalVar.mainwindow = this;
            this.KeyDown += new KeyEventHandler(tb_KeyDown);
            //TextCompositionManager.AddTextInputHandler(this, new TextCompositionEventHandler(OnTextComposition));
            //createGlobalCube_drawone();
            //return;
            if (GlobalVar.USEVERSION2 == true)
            {
                ExecuteFromMainWindow_v2();
            }
            else
            {
                ExecuteFromMainWindow();
            }
            Mycamera = this.camera;
            //DrawBar zz = new DrawBar(1,new Point3D(0, 0, 0), new Point3D(0, 10, 0),2);
            //AddArrow(10, 10);
            Vector3D vectorOfLook = new Vector3D();
            Vector3D updirect = Mycamera.UpDirection; ;
            /// remember to uncomment to reinstate axis and camera box
            Point3D pointOnXPLane = DrawLines(ref vectorOfLook);
            //AddCube(0, 0, pointOnXPLane, vectorOfLook);
            

        }
        private void testRequiredOrbitalParameters ()
        {
            double aaphelion = 152.10 * Math.Pow(10, 6);
            double aperihelion = 147.10 * Math.Pow(10, 6);
            double atilt =23.4;
            double aeccentricity = 0.01671123;
            double ainclination=7;
            double alongitudeoftheascendingnode = 348.73936;
            double aargumentofperiapsis = 114.20783;
            double aasemimajoraxis = 149.6 * Math.Pow(10, 6);
            string areferencePlane="";
            RequiredOrbitalParameters RO = new RequiredOrbitalParameters(aasemimajoraxis,aaphelion, aperihelion, atilt, aeccentricity, ainclination, alongitudeoftheascendingnode, aargumentofperiapsis, areferencePlane);
            EarthOrbit EO=RO.calculateOrbitRelatedVectors(Planets["sun"], new Vector3D(-1, 0, 0));
            EarthOrbit ggg = EO;

        }

        private void testPlanes ()
        {
            Vector3D n = new Vector3D(2, 4, 2);
            Plane pl = new Plane(n, new Point3D(1,3,1));
            EquationParameter EP=new EquationParameter();
            double res;
            EP.X=5.0;
            EP.Z=8.0;
            pl.Equation.solveForRemainingVariable(EP,out res);
            Vector3D[] basis = BASIS.findOrthonormalBasisToVector(n);
            object r = basis;
        }
        private void testGramSchmidt()
        {
            Vector3D x1 = new Vector3D(3,6,0);
            Vector3D x2 = new Vector3D(1,2,2);
            GramSchmidt GS = new GramSchmidt(false,x1,x2);
            object r = GS;





        }
        private Point3D DrawLines(ref Vector3D vectorOfLook)
        {
            Vector3D CamposOrig = (Vector3D)Mycamera.Position;

            Vector3D lookdirection = Mycamera.LookDirection;
            lookdirection.Normalize();
            Vector3D Campos = CamposOrig * .75;
            AddVisualAxix(new Point3D(-1000, 0, 0), new Point3D(1000, 0, 0), System.Windows.Media.Brushes.Green);
            AddVisualAxix(new Point3D(0, -1000, 0), new Point3D(0, 1000, 0), System.Windows.Media.Brushes.Red);
            //AddVisualAxix(new Point3D(0, -150, 0), new Point3D(0, 150, 0), System.Windows.Media.Brushes.SteelBlue);
            AddVisualAxix(new Point3D(0, 0, -1000), new Point3D(0, 0, 1000), System.Windows.Media.Brushes.White);
            Vector3D distnace;
            Vector3D projection;
            Vector3D normal= new Vector3D(0, 0, 1);
            projection = VectorProjections.calculateProjectionOfVectorIntoPlane(lookdirection, normal, out distnace);
            double realdistance = CamposOrig.Z;
            Vector3D realProjectionAlt = (realdistance / distnace.Length) * projection;
            //Vector3D realprojection = projection * realdistance;
            //Vector3D fianlVector = realprojection + distnace * realdistance;
            Vector3D  fianlVector = -lookdirection * (realdistance / distnace.Length);
            Point3D pointOnXPLane = (Point3D)(CamposOrig - fianlVector);
            AddVisualAxix((Point3D)CamposOrig-0.1*fianlVector, pointOnXPLane , System.Windows.Media.Brushes.Yellow);
            //AddVisualAxix((Point3D)CamposOrig - 0.0005 * fianlVector, new Point3D(0, 0, 0), System.Windows.Media.Brushes.Pink);
            //AddVisualAxix((Point3D)CamposOrig - 0.0005 * fianlVector, (Point3D)CamposOrig - 0.1 * fianlVector, System.Windows.Media.Brushes.Tan);

            vectorOfLook = fianlVector;
            return pointOnXPLane;
        }
        private void AddVisualAxix(Point3D A, Point3D B, System.Windows.Media.Brush color)
        {


            //MyModel.Transform = myRotateTransform;


            Point3D Campos = Mycamera.Position;
            Vector3D lookdirection = Mycamera.LookDirection / Mycamera.LookDirection.Length;

            double posscale = 10;
            //Matrix3D ofset = RotationMatrix.ReturnMatrix3DToOffsetByXYandZ_static(Campos.X / posscale, Campos.Y / posscale, Campos.Z / posscale);
            //ofset = RotationMatrix.ReturnMatrix3DToOffsetByXYandZ_static(0, 500, 700);
            //horizontal.Append(vertical);
            //MatrixTransform3D MT3D = new MatrixTransform3D(horizontal);

            DrawBar arrow = new DrawBar(2, A, B, 10, color);
            Model3DGroup  linegroup = arrow.group;
            Transform3DGroup transformGroupArrow = new Transform3DGroup();
            Transform3DGroup transformGroupArrow2 = new Transform3DGroup();
            double scalemeasure = 1;
            Matrix3D ScaleMArix = new Matrix3D(scalemeasure, 0, 0, 0,
                                              0, scalemeasure, 0, 0,
                                              0, 0, scalemeasure, 0,
                                             0, 0, 0, 1);
            ScaleTransform3D s = new ScaleTransform3D(scalemeasure, scalemeasure, scalemeasure);
            //AxisAngleRotation3D ax3d = new AxisAngleRotation3D(new Vector3D(0, 1, 0), angletoHorizontal + 0);
            //AxisAngleRotation3D ax3dv = new AxisAngleRotation3D(new Vector3D(1, 0, 0), angletoVertical + 0);
            //RotateTransform3D myRotateTransform = new RotateTransform3D(ax3d);
            //RotateTransform3D myRotateTransformv = new RotateTransform3D(ax3dv);
            TranslateTransform3D t = new TranslateTransform3D(Campos.X / posscale, Campos.Y / posscale, Campos.Z / posscale);
            //ScaleMArix.Append(horizontal);
            //ScaleMArix.Append(vertical);
            //ScaleMArix.Append(ofset);
            transformGroupArrow.Children.Add(new MatrixTransform3D(ScaleMArix));
            transformGroupArrow2.Children.Add(s);

            //transformGroupArrow2.Children.Add(myRotateTransform);
            //transformGroupArrow2.Children.Add(myRotateTransformv);
            //transformGroupArrow2.Children.Add(t);
            linegroup.Transform = transformGroupArrow;
            //arrowbodyandhead.Transform = transformGroupArrow2;
            viewport.Children.Add(new ModelVisual3D() { Content = linegroup });

        }
        private Model3DGroup ReturnModel3DGroupPerLine(double diameter, Point3D A, Point3D B, System.Windows.Media.Brush color)
        {


            //MyModel.Transform = myRotateTransform;


            //Point3D Campos = Mycamera.Position;
            //Vector3D lookdirection = Mycamera.LookDirection / Mycamera.LookDirection.Length;

            //double posscale = 10;
            //Matrix3D ofset = RotationMatrix.ReturnMatrix3DToOffsetByXYandZ_static(Campos.X / posscale, Campos.Y / posscale, Campos.Z / posscale);
            //ofset = RotationMatrix.ReturnMatrix3DToOffsetByXYandZ_static(0, 500, 700);
            //horizontal.Append(vertical);
            //MatrixTransform3D MT3D = new MatrixTransform3D(horizontal);

            DrawBar arrow = new DrawBar(2, A, B, diameter, color);
            Model3DGroup linegroup = arrow.group;
            Transform3DGroup transformGroupArrow = new Transform3DGroup();
            Transform3DGroup transformGroupArrow2 = new Transform3DGroup();
            double scalemeasure = 1;
            Matrix3D ScaleMArix = new Matrix3D(scalemeasure, 0, 0, 0,
                                              0, scalemeasure, 0, 0,
                                              0, 0, scalemeasure, 0,
                                             0, 0, 0, 1);
            //ScaleTransform3D s = new ScaleTransform3D(scalemeasure, scalemeasure, scalemeasure);
            //AxisAngleRotation3D ax3d = new AxisAngleRotation3D(new Vector3D(0, 1, 0), angletoHorizontal + 0);
            //AxisAngleRotation3D ax3dv = new AxisAngleRotation3D(new Vector3D(1, 0, 0), angletoVertical + 0);
            //RotateTransform3D myRotateTransform = new RotateTransform3D(ax3d);
            //RotateTransform3D myRotateTransformv = new RotateTransform3D(ax3dv);
            //TranslateTransform3D t = new TranslateTransform3D(Campos.X / posscale, Campos.Y / posscale, Campos.Z / posscale);
            //ScaleMArix.Append(horizontal);
            //ScaleMArix.Append(vertical);
            //ScaleMArix.Append(ofset);
            transformGroupArrow.Children.Add(new MatrixTransform3D(ScaleMArix));
            //transformGroupArrow2.Children.Add(s);

            //transformGroupArrow2.Children.Add(myRotateTransform);
            //transformGroupArrow2.Children.Add(myRotateTransformv);
            //transformGroupArrow2.Children.Add(t);
            linegroup.Transform = transformGroupArrow;
            //arrowbodyandhead.Transform = transformGroupArrow2;
            //viewport.Children.Add(new ModelVisual3D() { Content = linegroup });
            return linegroup;

        }
        private void AddArrow(double addx, double addy)
        {


            //MyModel.Transform = myRotateTransform;


            Point3D Campos = Mycamera.Position;
            Vector3D lookdirection = Mycamera.LookDirection / Mycamera.LookDirection.Length;
            double angletoHorizontal = addx + (Math.Atan(lookdirection.X / lookdirection.Z)) / Math.PI * 180;
            double angletoVertical = addy + (Math.Atan(lookdirection.Z / lookdirection.Y)) / Math.PI * 180; ;
            Matrix3D horizontal = RotationMatrix.ReturnRotationMatrixAroundAxixComingOutOfOrigin_inactionstatic(new Vector3D(0, 1, 0), angletoHorizontal);
            Matrix3D vertical = RotationMatrix.ReturnRotationMatrixAroundAxixComingOutOfOrigin_inactionstatic(new Vector3D(1, 0, 0), angletoVertical);
            
            
            
            double posscale = 10;
            Matrix3D ofset = RotationMatrix.ReturnMatrix3DToOffsetByXYandZ_static(Campos.X / posscale, Campos.Y / posscale, Campos.Z / posscale);
             //ofset = RotationMatrix.ReturnMatrix3DToOffsetByXYandZ_static(0, 500, 700);
            //horizontal.Append(vertical);
            //MatrixTransform3D MT3D = new MatrixTransform3D(horizontal);
            
            Arrow3DShape arrow = new Arrow3DShape();
            arrowbodyandhead = arrow.group;
            Transform3DGroup transformGroupArrow = new Transform3DGroup();
            Transform3DGroup transformGroupArrow2 = new Transform3DGroup();
            double scalemeasure = 20;
            Matrix3D ScaleMArix = new Matrix3D(scalemeasure, 0, 0, 0,
                                              0, scalemeasure, 0, 0,
                                              0, 0, scalemeasure, 0,
                                             0, 0, 0, 1);
            ScaleTransform3D s = new ScaleTransform3D(scalemeasure, scalemeasure, scalemeasure);
            AxisAngleRotation3D ax3d = new AxisAngleRotation3D(new Vector3D(0, 1, 0), angletoHorizontal + 0);
            AxisAngleRotation3D ax3dv = new AxisAngleRotation3D(new Vector3D(1, 0, 0), angletoVertical + 0);
            RotateTransform3D myRotateTransform = new RotateTransform3D(ax3d);
            RotateTransform3D myRotateTransformv = new RotateTransform3D(ax3dv);
            TranslateTransform3D t = new TranslateTransform3D(Campos.X / posscale, Campos.Y / posscale, Campos.Z / posscale);
            ScaleMArix.Append(horizontal);
            ScaleMArix.Append(vertical);
            //ScaleMArix.Append(ofset);
            transformGroupArrow.Children.Add(new MatrixTransform3D(ScaleMArix));
            transformGroupArrow2.Children.Add(s);
            
            transformGroupArrow2.Children.Add(myRotateTransform);
            transformGroupArrow2.Children.Add(myRotateTransformv);
            //transformGroupArrow2.Children.Add(t);
            arrowbodyandhead.Transform = transformGroupArrow;
            //arrowbodyandhead.Transform = transformGroupArrow2;
            viewport.Children.Add(new ModelVisual3D() { Content = arrowbodyandhead }); 

        }
        private void AddCube(double addx, double addy, Point3D pointOnXPLane, Vector3D vectorOfLook)
        {


            //MyModel.Transform = myRotateTransform;


            Point3D Campos = Mycamera.Position;
            Vector3D lookdirection = Mycamera.LookDirection / Mycamera.LookDirection.Length;
            double angletoHorizontal = (addx + (Math.Atan(lookdirection.X / lookdirection.Z)) / Math.PI * 180);
            double baseSideSquared=Math.Pow(lookdirection.Z, 2) + Math.Pow(lookdirection.X, 2);
            double baselength = Math.Sqrt(baseSideSquared);
            double angletoVertical =-( addy + (Math.Atan(baselength / lookdirection.Y)) / Math.PI * 180); ;
            //angletoHorizontal = angletoVertical = 0;
            Matrix3D horizontal = RotationMatrix.ReturnRotationMatrixAroundAxixComingOutOfOrigin_inactionstatic(new Vector3D(0, 1, 0), angletoHorizontal);
            Matrix3D vertical = RotationMatrix.ReturnRotationMatrixAroundAxixComingOutOfOrigin_inactionstatic(new Vector3D(1, 0, 0), angletoVertical);
            


            double posscale =1.25;
            Matrix3D ofset = RotationMatrix.ReturnMatrix3DToOffsetByXYandZ_static(Campos.X / posscale+500, Campos.Y / posscale, Campos.Z / posscale);

            //ofset = RotationMatrix.ReturnMatrix3DToOffsetByXYandZ_static(0, 500, 700);
            //horizontal.Append(vertical);
            //MatrixTransform3D MT3D = new MatrixTransform3D(horizontal);

            Cube cube = new Cube();
            CubeFrame = cube.group;
            Transform3DGroup transformGroupArrow = new Transform3DGroup();
            Transform3DGroup transformGroupArrow2 = new Transform3DGroup();
            Transform3DGroup transformGroupArrow_forgettingMatrix3d = new Transform3DGroup();
            double scalemeasure = 100;
            Matrix3D ScaleMArix = new Matrix3D(scalemeasure, 0, 0, 0,
                                              0, scalemeasure, 0, 0,
                                              0, 0, scalemeasure, 0,
                                             0, 0, 0, 1);
            ScaleTransform3D s = new ScaleTransform3D(scalemeasure, scalemeasure, scalemeasure);
            AxisAngleRotation3D ax3d = new AxisAngleRotation3D(new Vector3D(0, 1, 0), angletoHorizontal + 0);
            
            transformGroupArrow_forgettingMatrix3d.Children.Add(new RotateTransform3D(ax3d));
            Matrix3D MToRotateXAxis = transformGroupArrow_forgettingMatrix3d.Value;
            Point3D p = new Point3D(1, 0, 0);
            
            Point3D xAxixRotated_point = MToRotateXAxis.Transform(p);
            Vector3D xAxixRotated = (Vector3D)(xAxixRotated_point);
                

            AxisAngleRotation3D ax3dv = new AxisAngleRotation3D(xAxixRotated, angletoVertical + 0);
            RotateTransform3D myRotateTransform = new RotateTransform3D(ax3d);
            RotateTransform3D myRotateTransformv = new RotateTransform3D(ax3dv);
            TranslateTransform3D t0 = new TranslateTransform3D(-0.5, -0.5 , -0.5 );
            //vectorOfLook
            //TranslateTransform3D t = new TranslateTransform3D(Campos.X / posscale, Campos.Y / posscale, Campos.Z / posscale);
            TranslateTransform3D t = new TranslateTransform3D( pointOnXPLane.X ,  + pointOnXPLane.Y ,  + pointOnXPLane.Z);
            double dist=0.9;
            TranslateTransform3D t1 = new TranslateTransform3D(vectorOfLook.X*dist, +vectorOfLook.Y*dist, +vectorOfLook.Z*dist);
            ScaleMArix.Append(horizontal);
            ScaleMArix.Append(vertical);
            ScaleMArix.Append(ofset);
            transformGroupArrow.Children.Add(new MatrixTransform3D(ScaleMArix));
            transformGroupArrow2.Children.Add(t0);
            transformGroupArrow2.Children.Add(s);

            transformGroupArrow2.Children.Add(myRotateTransform);
            transformGroupArrow2.Children.Add(myRotateTransformv);
            double angleBetwenupdirectionPerprAndCurrent = findAngleBetweenTwoVectors(Mycamera.UpDirection, new Vector3D(0, 1, 0)); ;
            AxisAngleRotation3D updirectionrelated = new AxisAngleRotation3D(vectorOfLook, -angleBetwenupdirectionPerprAndCurrent/Math.PI*180);
            RotateTransform3D updirectionrelated_otateTransform3D = new RotateTransform3D(updirectionrelated);
            transformGroupArrow2.Children.Add(updirectionrelated_otateTransform3D);

            transformGroupArrow2.Children.Add(t);
            transformGroupArrow2.Children.Add(t1);
            CubeFrame.Transform = transformGroupArrow;
            CubeFrame.Transform = transformGroupArrow2;
            viewport.Children.Add(new ModelVisual3D() { Content = CubeFrame });

        }
        private static double findAngleBetweenTwoVectors (Vector3D v1,Vector3D v2      )
        {
            Vector3D v1n = v1;
            Vector3D v2n = v2;
            v1n.Normalize();
            v2n.Normalize();
            double cosAlpha = Vector3D.DotProduct(v1n,v2n);
            return Math.Acos(cosAlpha);
        }
        public static void createGlobalCube ()
        {


            //MyModel.Transform = myRotateTransform;



            //ofset = RotationMatrix.ReturnMatrix3DToOffsetByXYandZ_static(0, 500, 700);
            //horizontal.Append(vertical);
            //MatrixTransform3D MT3D = new MatrixTransform3D(horizontal);
            Transform3DGroup transformGroupArrow = new Transform3DGroup();
            Cube cube = new Cube();
            globalTrialCube = cube.group;


            double scalemeasure = 100;
            Matrix3D ScaleMArix = new Matrix3D(scalemeasure, 0, 0, 0,
                                             0, scalemeasure, 0, 0,
                                             0, 0, scalemeasure, 0,
                                            0, 0, 0, 1);
            //ScaleTransform3D s = new ScaleTransform3D(scalemeasure, scalemeasure, scalemeasure);
            ScaleTransform3D s = new ScaleTransform3D(1,1, 1);

            transformGroupArrow.Children.Add(s);

            globalTrialCube.Transform = transformGroupArrow;
            globalMesh = cube.Mesh;
            globalviewport.Children.Add(new ModelVisual3D() { Content = globalTrialCube });
 
        }
        public  void createGlobalCube_drawone()
        {

            InitializeComponent();
            //MyModel.Transform = myRotateTransform;



            //ofset = RotationMatrix.ReturnMatrix3DToOffsetByXYandZ_static(0, 500, 700);
            //horizontal.Append(vertical);
            //MatrixTransform3D MT3D = new MatrixTransform3D(horizontal);
            Transform3DGroup transformGroupArrow = new Transform3DGroup();
            Cube cube = new Cube(true);
            CubeFrame = cube.group;


            double scalemeasure = 100;
            Matrix3D ScaleMArix = new Matrix3D(scalemeasure, 0, 0, 0,
                                             0, scalemeasure, 0, 0,
                                             0, 0, scalemeasure, 0,
                                            0, 0, 0, 1);
            //ScaleTransform3D s = new ScaleTransform3D(scalemeasure, scalemeasure, scalemeasure);
            ScaleTransform3D s = new ScaleTransform3D(1, 0.5, 1);

            transformGroupArrow.Children.Add(s);

            CubeFrame.Transform = transformGroupArrow;
            globalMesh = cube.Mesh;
            viewport.Children.Add(new ModelVisual3D() { Content = CubeFrame });

        }
        private void OnTextComposition(object sender, TextCompositionEventArgs e)
        {
            string key = e.Text;
            Console.WriteLine(e.Text);
            ConsoleManager.Show();

        }
        static void tb_KeyDown(object sender, KeyEventArgs e)
        {
            Console.WriteLine(e.Key);
            ConsoleManager.Show();
            Vector3D lookdirection = Mycamera.LookDirection / Mycamera.LookDirection.Length;
            double step = 0.5;
            double angleBetweenZAxisAndLookDirection_cos=0;
            double angleBetweenZAxisAndLookDirection_sin = 0;
            Vector3D ZDirctioon=new Vector3D(0, 0, 1);
            Vector3D connectingLokDirectionAndZAxis = new Vector3D();
            Vector3D DirectionOfSideMovement = new Vector3D();
            if (e.Key == Key.Up)
            {
                //MoveCameraAccordingly(lookdirection*step);
                //MoveEyeVerticallyWhileRamainingFizatedOnSamePoint(1);
                MoveEyeHorizontallyOrVerticallyWhileRamainingFizatedOnSamePoint(1, false);
            }
            else if (e.Key == Key.W)
            {
                MoveCameraAccordingly(lookdirection * step);
            }
            else if (e.Key == Key.Down)
            {
                //MoveCameraAccordingly(lookdirection * - step);
                //MoveEyeVerticallyWhileRamainingFizatedOnSamePoint(-1);
                MoveEyeHorizontallyOrVerticallyWhileRamainingFizatedOnSamePoint(-1, false);
            }
            else if (e.Key == Key.Z)
            {
                MoveCameraAccordingly(lookdirection * -step);
            }
            else if (e.Key == Key.I)
            {
                ChangeLookDirection(-1, false);
            }
            else if (e.Key == Key.J)
            {
                ChangeLookDirection(1, true);
            }
            else if (e.Key == Key.K)
            {
                ChangeLookDirection(-1, true);
            }
            else if (e.Key == Key.M)
            {
                ChangeLookDirection(1,false);
            }
            else if (e.Key == Key.Left)
            {
                //angleBetweenZAxisAndLookDirection_cos = Vector3D.DotProduct(lookdirection, ZDirctioon);
                //angleBetweenZAxisAndLookDirection_sin = Math.Pow(1 - Math.Pow(angleBetweenZAxisAndLookDirection_cos, 2), 0.5);
                //connectingLokDirectionAndZAxis = lookdirection.Length / angleBetweenZAxisAndLookDirection_cos * ZDirctioon - lookdirection;// this one is perpendicular as we built rigth angle triangle with the result as one of the perpendicular sides
                //connectingLokDirectionAndZAxis = lookdirection.Length * angleBetweenZAxisAndLookDirection_cos * ZDirctioon - lookdirection;
                //if (connectingLokDirectionAndZAxis.Length == 0) connectingLokDirectionAndZAxis = new Vector3D(0, 1, 0);
                //DirectionOfSideMovement = Vector3D.CrossProduct(connectingLokDirectionAndZAxis, lookdirection);
                //DirectionOfSideMovement.Normalize();
                //MoveCameraAccordingly(DirectionOfSideMovement);
                //MoveRightOrLeftOld(1);
                MoveEyeHorizontallyOrVerticallyWhileRamainingFizatedOnSamePoint(1, true);
            }
            else if (e.Key == Key.Right)
            {
                //angleBetweenZAxisAndLookDirection_cos = Vector3D.DotProduct(lookdirection, ZDirctioon);
                //angleBetweenZAxisAndLookDirection_sin = Math.Pow(1 - Math.Pow(angleBetweenZAxisAndLookDirection_cos, 2), 0.5);
                //connectingLokDirectionAndZAxis = lookdirection.Length / angleBetweenZAxisAndLookDirection_cos * ZDirctioon - lookdirection;// this one is perpendicular as we built rigth angle triangle with the result as one of the perpendicular sides
                //connectingLokDirectionAndZAxis = lookdirection.Length * angleBetweenZAxisAndLookDirection_cos * ZDirctioon - lookdirection;
                //if (connectingLokDirectionAndZAxis.Length == 0) connectingLokDirectionAndZAxis = new Vector3D(0, 1, 0);
                //DirectionOfSideMovement = Vector3D.CrossProduct(connectingLokDirectionAndZAxis, lookdirection);
                //DirectionOfSideMovement.Normalize();
                //MoveCameraAccordingly(-DirectionOfSideMovement);
                //MoveRightOrLeftOld(-1);
                MoveEyeHorizontallyOrVerticallyWhileRamainingFizatedOnSamePoint(-1, true);
            }
            else if (e.Key == Key.A)
            {
                //enter key is down
                spinXDegreesClockWise_v2(1);
            }
            else if (e.Key == Key.S)
            {
                spinXDegreesClockWise_v2(-1);
                //enter key is down
            }
            else if (e.Key == Key.D)
            {
                //enter key is down
            }
            else if (e.Key == Key.U)
            {
                //enter key is down
            }
            else if (e.Key == Key.L)
            {
                //enter key is down
            }
            else if (e.Key == Key.R)
            {
                //enter key is down
                if (displayRay == true)
                {
                    displayRay = false;
                }
                else
                {
                    displayRay = true;
                }
            }
        }
        private static void ChangeLookDirection(double XMovement_InDegress,bool horizontalMovementFlag)
        {
            //XMovement_InDegress = 5;
            angleaccumalated_delet = angleaccumalated_delet + XMovement_InDegress;
            double stam = 0;
            if (Math.Abs(angleaccumalated_delet) == 89 || Math.Abs(angleaccumalated_delet) == 90)
            {
                stam = stam + 1;
            }
            stam = stam + 1;
            Vector3D currentUpDirection = Mycamera.UpDirection;
            Vector3D currentLookDirection = Mycamera.LookDirection;
            Point3D currentPosition = Mycamera.Position;


            Vector3D distVector;
            Vector3D normVectorZXPlane = new Vector3D(0, 1, 0);
            Vector3D normVectorYXPlane = new Vector3D(0, 0, 1);
            Vector3D normVectorYZPlane = new Vector3D(1, 0, 0);
            Point3D pointLookedAtConstiutingRottaionCenter;
            bool res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYXPlane, out pointLookedAtConstiutingRottaionCenter);
            if (res == false)
            {
                res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYZPlane, out pointLookedAtConstiutingRottaionCenter);
                //normVectorZXPlane = normVectorYZPlane;
                if (res == false)
                {
                    res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorZXPlane, out pointLookedAtConstiutingRottaionCenter);
                    //normVectorZXPlane = normVectorYZPlane;
                }
            
            }
            Vector3D lookVectorTowardPlaneInQuestion = currentLookDirection;// -(currentPosition - pointLookedAtConstiutingRottaionCenter);
            Vector3D lookdirectionUntivector = lookVectorTowardPlaneInQuestion;
            lookdirectionUntivector.Normalize();
            //Vector3D[] basisOfPlanePerpendicularTolookDirectiom = BASIS.findOrthonormalBasisToVector(lookdirectionUntivector);
            Vector3D ProjOntoPlanePerpendicularToLookDirection_distVector;
            // ProjOntoPlanePerpendicularToLookDirection is the up direction of the camera -- if we move vertically we move in that direction, depending on the angle with respect to the center found above
            Vector3D ProjOntoPlanePerpendicularToLookDirection_original;
             ProjOntoPlanePerpendicularToLookDirection_original = VectorProjections.calculateProjectionOfVectorIntoPlane(currentUpDirection, lookdirectionUntivector, out  ProjOntoPlanePerpendicularToLookDirection_distVector);//
             Vector3D ProjOntoPlanePerpendicularToLookDirection = ProjOntoPlanePerpendicularToLookDirection_original;
            ProjOntoPlanePerpendicularToLookDirection.Normalize();
            Vector3D ProjOntoPlanePerpendicularToLookDirection_horizontaltherto = Vector3D.CrossProduct(ProjOntoPlanePerpendicularToLookDirection, lookdirectionUntivector);
            double XMovement_InRadians = XMovement_InDegress / 180 * Math.PI;
            double CosX = Math.Cos(XMovement_InRadians);
            double SinX = Math.Sin(XMovement_InRadians);
            Vector3D NewLookVector_unitvector;
            Vector3D NewLookVector;
            Vector3D theotherperpendicular;
            Vector3D up;
            if (horizontalMovementFlag == true)
            {
                NewLookVector_unitvector = -lookdirectionUntivector * CosX + ProjOntoPlanePerpendicularToLookDirection_horizontaltherto * SinX;
                NewLookVector_unitvector = -NewLookVector_unitvector;
                theotherperpendicular = Vector3D.CrossProduct(-lookdirectionUntivector, ProjOntoPlanePerpendicularToLookDirection_horizontaltherto);
                up = ProjOntoPlanePerpendicularToLookDirection;
            }
            else
            {
                NewLookVector_unitvector =- lookdirectionUntivector * CosX + ProjOntoPlanePerpendicularToLookDirection * SinX;
                NewLookVector_unitvector = -NewLookVector_unitvector;
                theotherperpendicular = Vector3D.CrossProduct(-lookdirectionUntivector ,  ProjOntoPlanePerpendicularToLookDirection);
                theotherperpendicular.Normalize();
                up = -Vector3D.CrossProduct(theotherperpendicular, NewLookVector_unitvector);
            }
            
            
            //up.Normalize();
            NewLookVector = NewLookVector_unitvector;// *lookVectorTowardPlaneInQuestion.Length;
            //Point3D newpos;
            //newpos = pointLookedAtConstiutingRottaionCenter - NewLookVector;

            Mycamera.LookDirection = NewLookVector;
            //Mycamera.Position = newpos;
            
            Vector3D newUpDirection;
            //Mycamera.UpDirection = ProjOntoPlanePerpendicularToLookDirection;
            Mycamera.UpDirection = up;
            if (Mycamera.UpDirection.Length < 0.03)
            //if (findDistanceBetween2Vector(lookdirectionUntivector, ProjOntoPlanePerpendicularToLookDirection_original) < 0.05)
            {
                if (horizontalMovementFlag == true)
                {
                    newUpDirection = ProjOntoPlanePerpendicularToLookDirection_horizontaltherto;
                } else
                {
                    newUpDirection = ProjOntoPlanePerpendicularToLookDirection;
                }
                Mycamera.UpDirection = newUpDirection;
            }
            //Console.WriteLine("final position: " +Mycamera.Position);
            //Console.WriteLine("final length: " + Mycamera.LookDirection.Length + "\n" + "poition length: " + ((Vector3D)Mycamera.Position).Length + "\n" + "look length :" + Mycamera.LookDirection + "\n" + "position :" + Mycamera.Position + "\nup direction unti vector:" + ProjOntoPlanePerpendicularToLookDirection);
            Console.WriteLine(Mycamera.LookDirection + "\n" + "position :" + Mycamera.Position + "\nup direction unti vector:" + Mycamera.UpDirection);
            CompositionTarget_Rendering_V2_inaction();
            

        }
        private static void spinXDegreesClockWise_v2(double XMovement_InDegress)
        {
            //XMovement_InDegress = 5;
            angleaccumalated_delet = angleaccumalated_delet + XMovement_InDegress;
            double stam = 0;
            if (Math.Abs(angleaccumalated_delet) == 89 || Math.Abs(angleaccumalated_delet) == 90)
            {
                stam = stam + 1;
            }
            stam = stam + 1;
            Vector3D currentUpDirection = Mycamera.UpDirection;
            Vector3D currentLookDirection = Mycamera.LookDirection;
            Point3D currentPosition = Mycamera.Position;


            Vector3D distVector;
            Vector3D normVectorZXPlane = new Vector3D(0, 1, 0);
            Vector3D normVectorYXPlane = new Vector3D(0, 0, 1);
            Vector3D normVectorYZPlane = new Vector3D(1, 0, 0);
            Point3D pointLookedAtConstiutingRottaionCenter;
            bool res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYXPlane, out pointLookedAtConstiutingRottaionCenter);
            if (res == false)
            {
                res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYZPlane, out pointLookedAtConstiutingRottaionCenter);
                //normVectorZXPlane = normVectorYZPlane;
                if (res == false)
                {
                    res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorZXPlane, out pointLookedAtConstiutingRottaionCenter);
                    //normVectorZXPlane = normVectorYZPlane;
                }

            }
            Vector3D lookVectorTowardPlaneInQuestion = -(currentPosition - pointLookedAtConstiutingRottaionCenter);
            Vector3D lookdirectionUntivector = lookVectorTowardPlaneInQuestion;
            lookdirectionUntivector.Normalize();
            //Vector3D[] basisOfPlanePerpendicularTolookDirectiom = BASIS.findOrthonormalBasisToVector(lookdirectionUntivector);


            Vector3D ProjOntoPlanePerpendicularToLookDirection_distVector;
            // ProjOntoPlanePerpendicularToLookDirection is the up direction of the camera -- if we move vertically we move in that direction, depending on the angle with respect to the center found above
            Vector3D ProjOntoPlanePerpendicularToLookDirection_original;
            ProjOntoPlanePerpendicularToLookDirection_original = VectorProjections.calculateProjectionOfVectorIntoPlane(currentUpDirection, lookdirectionUntivector, out  ProjOntoPlanePerpendicularToLookDirection_distVector);//
            Vector3D ProjOntoPlanePerpendicularToLookDirection = ProjOntoPlanePerpendicularToLookDirection_original;
            ProjOntoPlanePerpendicularToLookDirection.Normalize();
            Vector3D ProjOntoPlanePerpendicularToLookDirection_horizontaltherto = Vector3D.CrossProduct(ProjOntoPlanePerpendicularToLookDirection, lookdirectionUntivector);
            Vector3D currentUpDirection_unit = currentUpDirection;
            currentUpDirection_unit.Normalize();
            
            
            double XMovement_InRadians = XMovement_InDegress / 180 * Math.PI;
            double CosX = Math.Cos(XMovement_InRadians);
            double SinX = Math.Sin(XMovement_InRadians);
            Vector3D NewLookVector_unitvector;
            Vector3D NewLookVector;

            NewLookVector_unitvector = currentUpDirection_unit * CosX + ProjOntoPlanePerpendicularToLookDirection_horizontaltherto * SinX;
            Mycamera.UpDirection = NewLookVector_unitvector;

            //NewLookVector = NewLookVector_unitvector * lookVectorTowardPlaneInQuestion.Length;
            //Point3D newpos;
            //newpos = pointLookedAtConstiutingRottaionCenter - NewLookVector;

            //Mycamera.LookDirection = NewLookVector;
            //Mycamera.Position = newpos;
            //CompositionTarget_Rendering_V2_inaction();
            //Vector3D newUpDirection;
            //if (ProjOntoPlanePerpendicularToLookDirection_original.Length < 0.03)
            //{

            //    newUpDirection = ProjOntoPlanePerpendicularToLookDirection_horizontaltherto;

            //    Mycamera.UpDirection = newUpDirection;
            //}
            //Console.WriteLine("final position: " +Mycamera.Position);
            Console.WriteLine("final length: " + Mycamera.LookDirection.Length + "\n" + "poition length: " + ((Vector3D)Mycamera.Position).Length);

                


        }
        private static void MoveEyeHorizontallyOrVerticallyWhileRamainingFizatedOnSamePoint(double XMovement_InDegress,bool horizontalMovementFlag)
        {
            //XMovement_InDegress = 5;
            angleaccumalated_delet = angleaccumalated_delet + XMovement_InDegress;
            double stam = 0;
            if (Math.Abs(angleaccumalated_delet) == 89 || Math.Abs(angleaccumalated_delet) == 90)
            {
                stam = stam + 1;
            }
            stam = stam + 1;
            Vector3D currentUpDirection = Mycamera.UpDirection;
            Vector3D currentLookDirection = Mycamera.LookDirection;
            Point3D currentPosition = Mycamera.Position;


            Vector3D distVector;
            Vector3D normVectorZXPlane = new Vector3D(0, 1, 0);
            Vector3D normVectorYXPlane = new Vector3D(0, 0, 1);
            Vector3D normVectorYZPlane = new Vector3D(1, 0, 0);
            Point3D pointLookedAtConstiutingRottaionCenter;
            bool res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYXPlane, out pointLookedAtConstiutingRottaionCenter);
            if (res == false)
            {
                res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYZPlane, out pointLookedAtConstiutingRottaionCenter);
                //normVectorZXPlane = normVectorYZPlane;
                if (res == false)
                {
                    res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorZXPlane, out pointLookedAtConstiutingRottaionCenter);
                    //normVectorZXPlane = normVectorYZPlane;
                }
            
            }
            Vector3D lookVectorTowardPlaneInQuestion = -(currentPosition - pointLookedAtConstiutingRottaionCenter);
            Vector3D lookdirectionUntivector = lookVectorTowardPlaneInQuestion;
            lookdirectionUntivector.Normalize();
            //Vector3D[] basisOfPlanePerpendicularTolookDirectiom = BASIS.findOrthonormalBasisToVector(lookdirectionUntivector);
            Vector3D ProjOntoPlanePerpendicularToLookDirection_distVector;
            // ProjOntoPlanePerpendicularToLookDirection is the up direction of the camera -- if we move vertically we move in that direction, depending on the angle with respect to the center found above
            Vector3D ProjOntoPlanePerpendicularToLookDirection_original;
             ProjOntoPlanePerpendicularToLookDirection_original = VectorProjections.calculateProjectionOfVectorIntoPlane(currentUpDirection, lookdirectionUntivector, out  ProjOntoPlanePerpendicularToLookDirection_distVector);//
             Vector3D ProjOntoPlanePerpendicularToLookDirection = ProjOntoPlanePerpendicularToLookDirection_original;
            ProjOntoPlanePerpendicularToLookDirection.Normalize();
            Vector3D ProjOntoPlanePerpendicularToLookDirection_horizontaltherto = Vector3D.CrossProduct(ProjOntoPlanePerpendicularToLookDirection, lookdirectionUntivector);
            double XMovement_InRadians = XMovement_InDegress / 180 * Math.PI;
            double CosX = Math.Cos(XMovement_InRadians);
            double SinX = Math.Sin(XMovement_InRadians);
            Vector3D NewLookVector_unitvector;
            Vector3D NewLookVector;
            if (horizontalMovementFlag == true)
            {
                NewLookVector_unitvector = lookdirectionUntivector * CosX + ProjOntoPlanePerpendicularToLookDirection_horizontaltherto * SinX;
            }
            else
            {
                NewLookVector_unitvector = lookdirectionUntivector * CosX + ProjOntoPlanePerpendicularToLookDirection * SinX;
            }

            NewLookVector = NewLookVector_unitvector * lookVectorTowardPlaneInQuestion.Length;
            Point3D newpos;
            newpos = pointLookedAtConstiutingRottaionCenter - NewLookVector;

            Mycamera.LookDirection = NewLookVector;
            Mycamera.Position = newpos;
            CompositionTarget_Rendering_V2_inaction();
            Vector3D newUpDirection;
            Mycamera.UpDirection = ProjOntoPlanePerpendicularToLookDirection;
            if ( ProjOntoPlanePerpendicularToLookDirection_original.Length < 0.03)
            //if (findDistanceBetween2Vector(lookdirectionUntivector, ProjOntoPlanePerpendicularToLookDirection_original) < 0.05)
            {
                if (horizontalMovementFlag == true)
                {
                    newUpDirection = ProjOntoPlanePerpendicularToLookDirection_horizontaltherto;
                } else
                {
                    newUpDirection = ProjOntoPlanePerpendicularToLookDirection;
                }
                Mycamera.UpDirection = newUpDirection;
            }
            //Console.WriteLine("final position: " +Mycamera.Position);
            Console.WriteLine("final length: " + Mycamera.LookDirection.Length + "\n" + "poition length: " + ((Vector3D)Mycamera.Position).Length + "\n" + "look length :" + Mycamera.LookDirection + "\n" + "position :" + Mycamera.Position + "\nup direction unti vector:" + ProjOntoPlanePerpendicularToLookDirection);

            return;
            Vector3D Proj;
            if (findDistanceBetween2Vector(new Vector3D(1, 0, 0), currentUpDirection) < 0.05)
            {
                Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorYZPlane, out  distVector);//
            }
            else
            {
                Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorZXPlane, out  distVector);//
            }

            //Vector3D Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorYZPlane, out  distVector);//
            // if (distVector.Length == 0)
            // {
            //      Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorYZPlane, out  distVector);//
            // }
            //if (distVector.Length == 0)
            // {
            //     Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorYXPlane, out  distVector);//
            // } 
            Console.WriteLine(lookVectorTowardPlaneInQuestion);
            Point3D ProjPlanePoint = (Point3D)(-(Vector3D)(distVector - pointLookedAtConstiutingRottaionCenter));
            //Proj = Proj * lookVectorTowardPlaneInQuestion.Length;
            Vector3D projUnitVector = Proj;
            projUnitVector.Normalize();
            Vector3D secondLegeOfBasisSpin = Vector3D.CrossProduct(projUnitVector, currentUpDirection); //first one is Proj

            secondLegeOfBasisSpin.Normalize();


            Vector3D newVecotrDirectionOnXZPlane = projUnitVector * CosX + secondLegeOfBasisSpin * SinX;
            Vector3D newVecorComingOutOfXZPlane = Proj.Length * newVecotrDirectionOnXZPlane;
            Point3D newFinalPosition = ProjPlanePoint - newVecorComingOutOfXZPlane;
            Vector3D newLookPosition = -(newFinalPosition - pointLookedAtConstiutingRottaionCenter);

            Mycamera.LookDirection = newLookPosition;
            Mycamera.Position = newFinalPosition;
            //Console.WriteLine(newLookPosition);
            //Console.WriteLine(newFinalPosition);
            //Console.WriteLine(currentUpDirection);

            CompositionTarget_Rendering_V2_inaction();

        }


        private static Vector3D[] returnBasisOflookDirectionAndPerpendicularToItOnSamePlaneOfLookDirectionAndOriginalUpDirection(Vector3D lookDirectionOnToCenteredPoint,Vector3D updirection)
        {
            GramSchmidt GM = new GramSchmidt(true, lookDirectionOnToCenteredPoint, updirection);
            return new Vector3D[] { GM.OrthonormalBasis[0], GM.OrthonormalBasis[1] };
        }
        private static void MoveEyeVerticallyWhileRamainingFizatedOnSamePoint(double XMovement_InDegress)
        {
            //XMovement_InDegress = 90;
            angleaccumalated_delet = angleaccumalated_delet + XMovement_InDegress;
            double stam = 0;
            if (Math.Abs(angleaccumalated_delet) == 89 || Math.Abs(angleaccumalated_delet) == 90)
            {
                stam = stam + 1;
            }
            stam = stam + 1;
            Vector3D currentUpDirection = Mycamera.UpDirection;
            Vector3D currentLookDirection = Mycamera.LookDirection;
            Point3D currentPosition = Mycamera.Position;


            Vector3D distVector;
            Vector3D normVectorZXPlane = new Vector3D(0, 1, 0);
            Vector3D normVectorYXPlane = new Vector3D(0, 0, 1);
            Vector3D normVectorYZPlane = new Vector3D(1, 0, 0);
            Point3D pointLookedAtConstiutingRottaionCenter;
            bool res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYXPlane, out pointLookedAtConstiutingRottaionCenter);
            if (res == false)
            {
                res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYZPlane, out pointLookedAtConstiutingRottaionCenter);

            }
            Vector3D lookVectorTowardPlaneInQuestion = -(currentPosition - pointLookedAtConstiutingRottaionCenter);
            Vector3D Proj;
            if (findDistanceBetween2Vector(new Vector3D(1, 0, 0), currentUpDirection) < 0.05)
            {
                Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorYZPlane, out  distVector);//
            }
            else
            {
                Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorZXPlane, out  distVector);//
            }

            Vector3D[] basis= returnBasisOflookDirectionAndPerpendicularToItOnSamePlaneOfLookDirectionAndOriginalUpDirection(-lookVectorTowardPlaneInQuestion, currentUpDirection);

            Console.WriteLine(lookVectorTowardPlaneInQuestion);
            Point3D ProjPlanePoint = (Point3D)(-(Vector3D)(distVector - pointLookedAtConstiutingRottaionCenter));

            Vector3D projUnitVector = Proj;
            projUnitVector.Normalize();
            Vector3D secondLegeOfBasisSpin = Vector3D.CrossProduct(projUnitVector, currentUpDirection); //first one is Proj



            double XMovement_InRadians = XMovement_InDegress / 180 * Math.PI;
            double CosX = Math.Cos(XMovement_InRadians);
            double SinX = Math.Sin(XMovement_InRadians);
            Vector3D newVecotrDirectionOnXZPlane = basis[0] * CosX + basis[1] * SinX;
            Vector3D newVecorComingOutOfXZPlane = lookVectorTowardPlaneInQuestion.Length * newVecotrDirectionOnXZPlane;
            Point3D newFinalPosition = newVecorComingOutOfXZPlane - pointLookedAtConstiutingRottaionCenter;
            Vector3D newLookPosition = -(newFinalPosition - pointLookedAtConstiutingRottaionCenter);

            Mycamera.LookDirection = newLookPosition;
            Mycamera.Position = newFinalPosition;
            Vector3D lookVectorTowardPlaneInQuestion_normalized = lookVectorTowardPlaneInQuestion;
            lookVectorTowardPlaneInQuestion_normalized.Normalize();
            if (findDistanceBetween2Vector(-lookVectorTowardPlaneInQuestion_normalized, currentUpDirection) < 0.025)
            {
                Mycamera.UpDirection = basis[1];// findPerpendiculartoUpDirection();
                //if (XMovement_InDegress < 0) Mycamera.UpDirection = Mycamera.UpDirection - Mycamera.UpDirection;
                currentUpDirection = Mycamera.UpDirection;
            }
            if (findDistanceBetween2Vector(lookVectorTowardPlaneInQuestion_normalized, currentUpDirection) < 0.025)
            {
                Mycamera.UpDirection = basis[1];// findPerpendiculartoUpDirection();
                //if (XMovement_InDegress < 0) Mycamera.UpDirection = Mycamera.UpDirection - Mycamera.UpDirection;
            }
            CompositionTarget_Rendering_V2_inaction();

        }
        private static Vector3D findPerpendiculartoUpDirection() 
        {
            return new Vector3D();
        }
        private static void MoveRightOrLeftOld(double direction)
        {
            MoveEyeHorizontallyWhileRamainingFizatedOnSamePoint_V0(-1.0 * direction);
            return;
            Vector3D lookdirection = Mycamera.LookDirection / Mycamera.LookDirection.Length;
   
            double angleBetweenZAxisAndLookDirection_cos = 0;
            double angleBetweenZAxisAndLookDirection_sin = 0;
            Vector3D ZDirctioon = new Vector3D(0, 0, 1);
            Vector3D connectingLokDirectionAndZAxis = new Vector3D();
            Vector3D DirectionOfSideMovement = new Vector3D();
            angleBetweenZAxisAndLookDirection_cos = Vector3D.DotProduct(lookdirection, ZDirctioon);
            angleBetweenZAxisAndLookDirection_sin = Math.Pow(1 - Math.Pow(angleBetweenZAxisAndLookDirection_cos, 2), 0.5);
            connectingLokDirectionAndZAxis = lookdirection.Length / angleBetweenZAxisAndLookDirection_cos * ZDirctioon - lookdirection;// this one is perpendicular as we built rigth angle triangle with the result as one of the perpendicular sides
            connectingLokDirectionAndZAxis = lookdirection.Length * angleBetweenZAxisAndLookDirection_cos * ZDirctioon - lookdirection;
            if (connectingLokDirectionAndZAxis.Length == 0) connectingLokDirectionAndZAxis = new Vector3D(0, 1, 0);
            DirectionOfSideMovement = Vector3D.CrossProduct(connectingLokDirectionAndZAxis, lookdirection);
            DirectionOfSideMovement.Normalize();
            MoveCameraAccordingly(DirectionOfSideMovement*direction);
        }
        private static void MoveEyeHorizontallyWhileRamainingFizatedOnSamePoint(double XMovement_InDegress)
        {
            angleaccumalated_delet = angleaccumalated_delet + XMovement_InDegress;
            double stam=0;
            if (Math.Abs(angleaccumalated_delet) == 89 || Math.Abs(angleaccumalated_delet) == 90)
            {
                stam = stam + 1;
            }
            stam = stam + 1;
            //XMovement_InDegress = 0;
            Vector3D currentUpDirection = Mycamera.UpDirection;
            Vector3D currentLookDirection = Mycamera.LookDirection;
            Point3D currentPosition = Mycamera.Position;


            Vector3D distVector;
            Vector3D normVectorZXPlane = new Vector3D(0, 1, 0);
            Vector3D normVectorYXPlane = new Vector3D(0, 0, 1);
            Vector3D normVectorYZPlane = new Vector3D(1, 0, 0);
            Point3D pointLookedAtConstiutingRottaionCenter;
            bool res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYXPlane, out pointLookedAtConstiutingRottaionCenter);
            if (res == false)
            {
                res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYZPlane, out pointLookedAtConstiutingRottaionCenter);
            }
            Vector3D lookVectorTowardPlaneInQuestion = -(currentPosition - pointLookedAtConstiutingRottaionCenter);

            Vector3D Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorZXPlane, out  distVector);
            Point3D ProjPlanePoint = (Point3D)(-(Vector3D)(distVector - pointLookedAtConstiutingRottaionCenter));
            //Proj = Proj * lookVectorTowardPlaneInQuestion.Length;
            Vector3D projUnitVector = Proj;
            projUnitVector.Normalize();
            Vector3D firstBasisLeg_interim = new Vector3D(-currentUpDirection.Y, currentUpDirection.X,0);
            Vector3D aa=firstBasisLeg_interim;
            Vector3D secondLegeOfBasisSpin = Vector3D.CrossProduct(projUnitVector, normVectorZXPlane); //first one is Proj
            //secondLegeOfBasisSpin = firstBasisLeg_interim;


            double XMovement_InRadians = XMovement_InDegress / 180 * Math.PI;
            double CosX = Math.Cos(XMovement_InRadians);
            double SinX = Math.Sin(XMovement_InRadians);
            Vector3D newVecotrDirectionOnXZPlane = projUnitVector * CosX + secondLegeOfBasisSpin * SinX;
            Vector3D newVecorComingOutOfXZPlane = Proj.Length * newVecotrDirectionOnXZPlane;
            Point3D newFinalPosition = ProjPlanePoint - newVecorComingOutOfXZPlane;
            Vector3D newLookPosition = -(newFinalPosition - pointLookedAtConstiutingRottaionCenter);

            Mycamera.LookDirection = newLookPosition;
            Mycamera.Position = newFinalPosition;


        }
        private static void MoveEyeHorizontallyWhileRamainingFizatedOnSamePoint_V0(double XMovement_InDegress)
        {
            angleaccumalated_delet = angleaccumalated_delet + XMovement_InDegress;
            double stam = 0;
            if (Math.Abs(angleaccumalated_delet) == 89 || Math.Abs(angleaccumalated_delet) == 90)
            {
                stam = stam + 1;
            }
            stam = stam + 1;
            Vector3D currentUpDirection = Mycamera.UpDirection;
            Vector3D currentLookDirection = Mycamera.LookDirection;
            Point3D currentPosition = Mycamera.Position;


            Vector3D distVector;
            Vector3D normVectorZXPlane = new Vector3D(0, 1, 0);
            Vector3D normVectorYXPlane = new Vector3D(0, 0, 1);
            Vector3D normVectorYZPlane = new Vector3D(1, 0, 0);
            Point3D pointLookedAtConstiutingRottaionCenter;
            bool res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYXPlane, out pointLookedAtConstiutingRottaionCenter);
            if (res == false)
            {
                res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYZPlane, out pointLookedAtConstiutingRottaionCenter);
                //normVectorZXPlane = normVectorYZPlane;
            }
            Vector3D lookVectorTowardPlaneInQuestion = -(currentPosition - pointLookedAtConstiutingRottaionCenter);
            Vector3D Proj;
            if (findDistanceBetween2Vector(new Vector3D(1, 0, 0), currentUpDirection) < 0.05)
            {
                Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorYZPlane, out  distVector);//
            }
            else
            {
                Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorZXPlane, out  distVector);//
            }
            
            //Vector3D Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorYZPlane, out  distVector);//
           // if (distVector.Length == 0)
           // {
           //      Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorYZPlane, out  distVector);//
           // }
           //if (distVector.Length == 0)
           // {
           //     Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorYXPlane, out  distVector);//
           // } 
            Console.WriteLine(lookVectorTowardPlaneInQuestion);
            Point3D ProjPlanePoint = (Point3D)(-(Vector3D)(distVector - pointLookedAtConstiutingRottaionCenter));
            //Proj = Proj * lookVectorTowardPlaneInQuestion.Length;
            Vector3D projUnitVector = Proj;
            projUnitVector.Normalize();
            Vector3D secondLegeOfBasisSpin = Vector3D.CrossProduct(projUnitVector,currentUpDirection); //first one is Proj

            secondLegeOfBasisSpin.Normalize();

            double XMovement_InRadians = XMovement_InDegress / 180 * Math.PI;
            double CosX = Math.Cos(XMovement_InRadians);
            double SinX = Math.Sin(XMovement_InRadians);
            Vector3D newVecotrDirectionOnXZPlane = projUnitVector * CosX + secondLegeOfBasisSpin * SinX;
            Vector3D newVecorComingOutOfXZPlane = Proj.Length * newVecotrDirectionOnXZPlane;
            Point3D newFinalPosition = ProjPlanePoint - newVecorComingOutOfXZPlane;
            Vector3D newLookPosition = -(newFinalPosition - pointLookedAtConstiutingRottaionCenter);

            Mycamera.LookDirection = newLookPosition;
            Mycamera.Position = newFinalPosition;
            //Console.WriteLine(newLookPosition);
            //Console.WriteLine(newFinalPosition);
            //Console.WriteLine(currentUpDirection);

            CompositionTarget_Rendering_V2_inaction();

        }
        private static double findDistanceBetween2Vector(Vector3D a , Vector3D b)
        {
            double distSquare = Math.Pow(a.X-b.X,2)+Math.Pow(a.Y-b.Y,2)+Math.Pow(a.Z-b.Z,2);
            return Math.Sqrt(distSquare);
        }
        private static void MoveEyeHorizontallyWhileRamainingFizatedOnSamePoint_scrap(double XMovement_InDegress)
        {
            Vector3D currentUpDirection = Mycamera.UpDirection;
            Vector3D currentLookDirection = Mycamera.LookDirection;
            Point3D currentPosition= Mycamera.Position;
          

            Vector3D distVector;
            Vector3D normVectorZXPlane=new Vector3D(0,1,0);
            Point3D pointOnXZPlane;
            bool res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorZXPlane, out pointOnXZPlane);
            Vector3D lookVectorTowardPositionFromXZPalne = currentPosition - pointOnXZPlane;

            Vector3D Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPositionFromXZPalne, normVectorZXPlane, out  distVector);
            Proj = Proj * lookVectorTowardPositionFromXZPalne.Length;
            double angleBetweenXZPlaneAndlookVectorTowardPositionFromXZPalne = findAngleBetweenTwoVectors(lookVectorTowardPositionFromXZPalne, Proj);
            double CurrentAngle_indegrees = angleBetweenXZPlaneAndlookVectorTowardPositionFromXZPalne / Math.PI * 180;
            double CosCurrentAngleE = Math.Cos(angleBetweenXZPlaneAndlookVectorTowardPositionFromXZPalne);
            double SinCurrentAngleE = Math.Sin(angleBetweenXZPlaneAndlookVectorTowardPositionFromXZPalne);
            //Vector3D perpendicularToProjAndLookDirectionTowardsXZPalne = Vector3D.CrossProduct(Proj, normVectorZXPlane);
            Vector3D normVectorZXPlane_adjusted = Proj / Proj.Length * SinCurrentAngleE + normVectorZXPlane * CosCurrentAngleE;
            Vector3D perpendicularToInitialUpdirection_adjustedAndLookDirectionToXZPalnePoint = Vector3D.CrossProduct(normVectorZXPlane_adjusted, lookVectorTowardPositionFromXZPalne);



            double CurrentAngleBetweenUpdirectionAndPerPendicular = findAngleBetweenTwoVectors(currentUpDirection, normVectorZXPlane);
            double Cos = Math.Cos(CurrentAngleBetweenUpdirectionAndPerPendicular);
            double Sin = Math.Sin(CurrentAngleBetweenUpdirectionAndPerPendicular);
            
            
            double CurrentAngleBetweenUpdirectionAndPerPendicular_indegrees = CurrentAngleBetweenUpdirectionAndPerPendicular / Math.PI * 180;
            Vector3D Final_normVectorZXPlane_adjusted; // to adjust if camera is not really held tilted using normVectorZXPlane_adjusted and perpendicularToInitialUpdirection_adjustedAndLookDirectionToXZPalnePoint as basis
            perpendicularToInitialUpdirection_adjustedAndLookDirectionToXZPalnePoint.Normalize();
            normVectorZXPlane_adjusted.Normalize();

            // this is the RELATIVE up diretion to the look direction, allowing for the camera to tilt
            Final_normVectorZXPlane_adjusted = normVectorZXPlane_adjusted * Cos + perpendicularToInitialUpdirection_adjustedAndLookDirectionToXZPalnePoint * Sin;
            // to acomplish movemnt of lok directiob by x degress, we need to use an ortonormal basis, one of its vector is the lookVectorTowardPositionFromXZPalne 
            // and the other is perpendicular to Final_normVectorZXPlane_adjusted and look direction lookVectorTowardPositionFromXZPalne

            Vector3D perpOfFinalBasis = Vector3D.CrossProduct(lookVectorTowardPositionFromXZPalne, Final_normVectorZXPlane_adjusted);
            perpOfFinalBasis.Normalize();
            Vector3D normalized_lookVectorTowardPositionFromXZPalne = lookVectorTowardPositionFromXZPalne;
            normalized_lookVectorTowardPositionFromXZPalne.Normalize();
            perpOfFinalBasis.Normalize();
            double XMovement_InRadians = XMovement_InDegress/180*Math.PI;
            double CosX = Math.Cos(XMovement_InRadians);
            double SinX = Math.Sin(XMovement_InRadians);
            Vector3D newVecorComingOutOfXZPlane = lookVectorTowardPositionFromXZPalne.Length * (normalized_lookVectorTowardPositionFromXZPalne * CosX + perpOfFinalBasis*SinX);
            Point3D newFinalPosition = pointOnXZPlane + newVecorComingOutOfXZPlane;
        
        }
        private static void spinXDegreesClockWise(double x)
        {

            Vector3D currentUpDirection = Mycamera.UpDirection;
            Vector3D perpenVector = new Vector3D(0, 1, 0);
            double CurrentAngle = findAngleBetweenTwoVectors(currentUpDirection,perpenVector);
            double CurrentAngle_indegrees = CurrentAngle / Math.PI * 180;
            if (currentUpDirection.X<0) CurrentAngle_indegrees=CurrentAngle_indegrees*-1;
            double finalAngleInDegrees = CurrentAngle_indegrees + x;
            double finalAngleInRAdians = finalAngleInDegrees / 180 * Math.PI;
  
            double CosCurrentAngleE = Math.Cos(finalAngleInRAdians);
            double SinCurrentAngleE = Math.Sin(finalAngleInRAdians);
            Vector3D finalUpDirection = new Vector3D(SinCurrentAngleE, CosCurrentAngleE, 0);
            Mycamera.UpDirection = finalUpDirection;

        }
        protected bool FindVectorCreatingAngleAlphaWithGivenVector(double alpha, Vector3D VectorInquestion, double x, double y, double z, ref Vector3D ResultVector)
        {

            bool result = true;
            double dotproductpartial;

            if (!(VectorInquestion.Y == 0))// /(PlaneInWhichPointIsSought == "Y")
            {
                if (VectorInquestion.Y == 0)
                {
                    result = false;
                }
                else
                {
                    dotproductpartial = VectorInquestion.X * x + VectorInquestion.Z * z;
                    y = -dotproductpartial / VectorInquestion.Y;
                }


            }
            else if (!(VectorInquestion.X == 0))// (PlaneInWhichPointIsSought == "X")
            {
                if (VectorInquestion.X == 0)
                {
                    result = false;
                }
                else
                {
                    dotproductpartial = VectorInquestion.Y * y + VectorInquestion.Z * z;
                    x = -dotproductpartial / VectorInquestion.X;
                }


            }
            else if (!(VectorInquestion.Z == 0))
            {
                if (VectorInquestion.Z == 0)
                {
                    result = false;
                }
                else
                {

                    dotproductpartial = VectorInquestion.X * x + VectorInquestion.Y * y;
                    z = -dotproductpartial / VectorInquestion.Z;
                }

            }


            if (result == true) ResultVector = new Vector3D(x, y, z);
            return result;

        }
        private static void MoveCameraAccordingly(Vector3D delta)
        {
            Vector3D lookdirection = Mycamera.LookDirection;
            Mycamera.Position = Mycamera.Position + delta;
        }
        private void ExecuteFromMainWindow()
        {
            PopulatePlanets();
            testPlaneIntersection();
            SomeTestWithVectors();
            testRotationMatrix();
            InitializeComponent();
            Point3D pointOfExternalAxis = new Point3D(0, 0, 0);
            Vector3D vectorOfExternalAxis = new Vector3D(0, 1, 0);
            Build3DTube(ReturnSHapesList(pointOfExternalAxis), pointOfExternalAxis, vectorOfExternalAxis);
            //CompositionTarget.Rendering += CompositionTarget_Rendering;
            CompositionTarget.Rendering += CompositionTarget_Rendering_V2;
        }
        private void ExecuteFromMainWindow_v2()
        {
            PopulatePlanets();
            testPlaneIntersection();
            SomeTestWithVectors();
            testRotationMatrix();
            InitializeComponent();
            Point3D pointOfExternalAxis = new Point3D(0, 0, 0);
            Vector3D vectorOfExternalAxis = new Vector3D(0, 1, 0);
            Build3DTube_v2(ReturnSHapesList_v2(pointOfExternalAxis), pointOfExternalAxis, vectorOfExternalAxis);
            //CompositionTarget.Rendering += CompositionTarget_Rendering;
            CompositionTarget.Rendering += CompositionTarget_Rendering_V2;
            testRequiredOrbitalParameters();
        }
        private void PopulatePlanets()
        {
            // Planet Sun
            PlanetPropeties sun = new PlanetPropeties();

            //public const double VisualToRealRAtio=10.0;

            Point3D sunCenter = new Point3D(0, 0, 0);

            OrbitAroundData OAD_forsun = new OrbitAroundData();
            sun.mass = 1.9885*Math.Pow(10,30);
            sun.aphelion = 0;
            sun.perihelion = 0;
            sun.SpeedAtperihelion = 0;
            sun.SpeedAtaphelion = 0;
            sun.OrbitSemiMinorAxis = 0;
            sun.OrbitSemiMajorAxis = 0;
            sun.Tilt = 0;
            sun.eccentricity = 0;
            //sun.MyOrbitFocalPlanet = null;
            sun.PlaneIntersectingOrbit = new Plane(new Vector3D(0, 0, 1), sunCenter);
            sun.initilaLocation = sunCenter;
            sun.MyOrbitFocalPlanet = OAD_forsun;
            sun.mormalToPlanetEquator = new Vector3D(0, 1, 0);
            double sunMagnifier = 50;
            PlanetDimensions sunDimesnions=new PlanetDimensions();

            sunDimesnions.PolarRadius=sunMagnifier*6.96 * Math.Pow(10,5)  ;
            sunDimesnions.EquaterialSemiMajorLength=sunDimesnions.PolarRadius ;
            sunDimesnions.EquaterialSemiMinorLength=sunDimesnions.PolarRadius;

            sun.Dimensions=sunDimesnions;

            // Planet Earth
            PlanetPropeties OtherEarth = new PlanetPropeties();
            OrbittingPlane myOrbittinPlaneOther = new OrbittingPlane();
            myOrbittinPlaneOther.planet = sun;
            myOrbittinPlaneOther.inclinationWithRespectToNormal = 0;
            OtherEarth.mass = 5.97219 * Math.Pow(10, 24);
            OtherEarth.aphelion = 152.10 * Math.Pow(10, 6);
            OtherEarth.perihelion = 147.10 * Math.Pow(10, 6);
            OtherEarth.SpeedAtperihelion = 0;
            OtherEarth.SpeedAtaphelion = 0;
            OtherEarth.OrbitSemiMinorAxis = 0;
            OtherEarth.eccentricity = 0.01671123;
            OtherEarth.OrbitSemiMajorAxis = 149.60 * Math.Pow(10, 6);
            OtherEarth.calculateSemiMinor();
            OtherEarth.Tilt = 23.4;
            OtherEarth.angleWithRespectToInitialLocation = 0;

            OrbitAroundData OADOther = new OrbitAroundData();
            OADOther.planet = sun;
            OADOther.intersectionPlane = new Plane(new Vector3D(1, 0, 0), new Point3D(0, 0, 0));
            OADOther.usePlanetLocatonAsPointInIntersectingLine = true;
            OADOther.OrbittingPlaneAngleDeterminator = myOrbittinPlaneOther;
            OtherEarth.MyOrbitFocalPlanet = OADOther;
            OtherEarth.myCurrentLocation = OtherEarth.initilaLocation;
            //Earth.PlaneIntersectingOrbit = new Plane(new Vector3D(0, 0, 1), sunCenter);
            double EarthDimMagnifierOther = 500.0;
            PlanetDimensions EarthDimesnionsOther = new PlanetDimensions();

            EarthDimesnionsOther.PolarRadius = EarthDimMagnifierOther * 6356.89;
            EarthDimesnionsOther.EquaterialSemiMajorLength = EarthDimMagnifierOther * 12756.3 / 2;
            EarthDimesnionsOther.EquaterialSemiMinorLength = EarthDimesnionsOther.EquaterialSemiMajorLength;

            OtherEarth.Dimensions = EarthDimesnionsOther;

            OtherEarth.CalculateSemiMajorAxixProjectionOntoReferncePlane();
            OtherEarth.calculateNormalToEquatorAndSemiMAjorAxisDirection();
            OtherEarth.CalcualteImaginaryCenter();
   
            OtherEarth.calculateSemiMinor();
            double angleInradians_otherearth = OtherEarth.angleWithRespectToInitialLocation / 180 * Math.PI;
            double cos_otherearth = Math.Cos(angleInradians_otherearth);
            double sin_otherearth = Math.Sin(angleInradians_otherearth);
            OtherEarth.initilaLocation = OtherEarth.imaginaryCenter + OtherEarth.SemiMajorAxisDirection * OtherEarth.OrbitSemiMajorAxis * cos_otherearth + OtherEarth.SemiMinorAxisDirection * OtherEarth.OrbitSemiMinorAxis * sin_otherearth;
            OtherEarth.CalcualteOtherFocalPoint();


            // Planet Earth
            PlanetPropeties Earth = new PlanetPropeties();
            OrbittingPlane myOrbittinPlane = new OrbittingPlane();
            myOrbittinPlane.planet = sun;
            myOrbittinPlane.inclinationWithRespectToNormal = 7;
            Earth.mass = 5.97219 * Math.Pow(10, 24);
            Earth.aphelion = 152.10 * Math.Pow(10, 6);
            Earth.perihelion = 147.10 * Math.Pow(10, 6);
            Earth.SpeedAtperihelion = 0;
            Earth.SpeedAtaphelion = 0;
            Earth.OrbitSemiMinorAxis = 0;
            Earth.eccentricity = 0.01671123;
            Earth.OrbitSemiMajorAxis = 149.60 * Math.Pow(10, 6);
            Earth.calculateSemiMinor();
            Earth.Tilt = 23.4;
            Earth.angleWithRespectToInitialLocation = 90;

            OrbitAroundData OAD = new OrbitAroundData();
            OAD.planet = sun;
            OAD.intersectionPlane = new Plane(new Vector3D(1, 0, 0), new Point3D(0, 0, 0));
            OAD.usePlanetLocatonAsPointInIntersectingLine = true;
            OAD.OrbittingPlaneAngleDeterminator = myOrbittinPlane;
            Earth.MyOrbitFocalPlanet = OAD;
            //Earth.PlaneIntersectingOrbit = new Plane(new Vector3D(0, 0, 1), sunCenter);
            double EarthDimMagnifier = 500.0;
            PlanetDimensions EarthDimesnions = new PlanetDimensions();

            EarthDimesnions.PolarRadius = EarthDimMagnifier * 6356.89;
            EarthDimesnions.EquaterialSemiMajorLength = EarthDimMagnifier * 12756.3 / 2;
            EarthDimesnions.EquaterialSemiMinorLength = EarthDimesnions.EquaterialSemiMajorLength;

            Earth.Dimensions = EarthDimesnions;

            Earth.CalculateSemiMajorAxixProjectionOntoReferncePlane();
            Earth.calculateNormalToEquatorAndSemiMAjorAxisDirection();
            Earth.CalcualteImaginaryCenter();
            
            Earth.calculateSemiMinor();
            double angleInradians = Earth.angleWithRespectToInitialLocation / 180 * Math.PI;
            double cos = Math.Cos(angleInradians);
            double sin = Math.Sin(angleInradians);
            Earth.initilaLocation = Earth.imaginaryCenter + Earth.SemiMajorAxisDirection * Earth.OrbitSemiMajorAxis * cos + Earth.SemiMinorAxisDirection * Earth.OrbitSemiMinorAxis * sin;
            Earth.CalcualteOtherFocalPoint();
            Earth.myCurrentLocation = Earth.initilaLocation;



            Planets["sun"] = sun;
            Planets["earth"] = Earth;
            Planets["earth2"] = OtherEarth;
        }
        private void testPlaneIntersection ()
        {
            Plane A = new Plane(new Vector3D(1,0,0),new Point3D(8,0,0));
            Plane B = new Plane(new Vector3D(0, 1, 0), new Point3D(0, 5, 0));
            LineVector LV = Plane.ReturnIntersectionLine(A, B);
            object k = LV;
        }
        private void SomeTestWithVectors()
        {
            Vector3D A = new Vector3D(-1,2,4);
            Vector3D V = new Vector3D(Math.Pow(0.5, 0.5), 1, Math.Pow(0.5, 0.5));
            Vector3D C = Vector3D.CrossProduct(V, A);
            Vector3D D = C;
            double a = Vector3D.DotProduct(C,V);
            double b = a;
            Vector3D A1 = new Vector3D(3, -3, 1);
            Vector3D V1 = new Vector3D(4, 9, 2);
            Vector3D C1 = Vector3D.CrossProduct(A1, V1);
            double D1 = Vector3D.DotProduct(C1,A1);
            double E1= D1;
            Vector3D normalOfPlaneInWhichTiltVectorProducesTiltAngleWithYAxix = new Vector3D(0, 0, -1);
            RequiredOrbitalParameters.calculateVernalPoint(Planets["earth"],normalOfPlaneInWhichTiltVectorProducesTiltAngleWithYAxix
                /* Planets["sun"]*/
                );
        }
        private void testRotationMatrix()
        {
            RotationMatrix MyMAtrix = new RotationMatrix(new Vector3D(0,1,0),new Point3D(0,0,0),1.25);
            Point3D P = new Point3D(2, 0, 0);
            Matrix3D A = MyMAtrix.ReturnRotationMatrix(90);
            
           Point3D final= A.Transform(P);
            Point3D s =final;
           Point3D z = (MyMAtrix.ReturnRotationMatrix(90)).Transform(P);
           Point3D f = (MyMAtrix.ReturnRotationMatrix(225)).Transform(P);
           Point3D r = A.Transform(f);
           Point3D v = A.Transform(r);
           Point3D m =v;

        }
        void testmatrix ()
        {
            Vector3D vector1 = new Vector3D(20, 30, 40);
            Vector3D vector2= new Vector3D(-200, -200, 40);
            Matrix3D matrix1 = new Matrix3D(10, 10, 10, 0, 20, 20, 20, 0, 30, 30, 30, 0, 5, 10, 15, 1);
            Vector3D vectorResult = new Vector3D();
            vectorResult = Vector3D.Multiply(vector1,matrix1);

            Point3D inM = new Point3D(250, 250, 0);

            Matrix3D Screen2World;
            //setup matrix for screen to world
            Screen2World = Matrix3D.Identity;

            Screen2World.Translate(new Vector3D(-200, -200, 0));
            Point3D outM = Point3D.Multiply(inM, Screen2World);
            Vector3D torres = Vector3D.Multiply(vector2, Screen2World); ;
        }
        private static  void CompositionTarget_Rendering_V2_inaction()
        {

            foreach (TubeGeometry TG in ShapesList)
            {







                TG.ExternalAxis_v2.CalculateMatrix();
                //TG.ExternalAxis.ApplyMAtrixToDixtionaryOfPoints(TG.Important3DPoints);

                Matrix3D RotationAroundExternalAxis = TG.ExternalAxis_v2.RotationAroundArbAxix.RotationMatrix3D;

                ((Transform3DGroup)(TG.GeometriesAndMtrices.mGeometry.Transform)).Children[0] = new MatrixTransform3D(RotationAroundExternalAxis);
                if (TG.DisplayPole == true)
                {
                    ((Transform3DGroup)(TG.GeometriesAndMtrices.mGeometry3.Transform)).Children[0] = new MatrixTransform3D(RotationAroundExternalAxis);
                }

                if (TG.DisplayCommet == true && 1 == 0)
                {
                    Matrix3D TranslateMatrixInitial2 = returnMatrixPertransformation(0, TG.TheAxis, new Vector3D(0, 0, 0));//selfrotatingangle was removed
                    ((Transform3DGroup)(TG.GeometriesAndMtrices.mGeometry2.Transform)).Children[0] = new MatrixTransform3D(TranslateMatrixInitial2);
                }
                if (TG.DisplayTrailFlag == true)
                {
                    TG.Trail.AddPointIfApplicable(TG.Important3DPoints_V2["CENTER POINT"].Now);
                }

                //SpinMode = SpinMode + 1;
                //if (SpinMode == 4) SpinMode = 0;
                //break;
                //displayRayIfNecessary();
            }

        }
        private  void displayRayIfNecessary()
        {
           
            if (displayRay == false)
            {
                if (!(lookDirectionRay == null))
                {
                    viewport.Children.Remove(lookDirectionRay);
                    lookDirectionRay = null;
                }
                
            }
            else
            {
                if (!(lookDirectionRay == null))
                {
                    viewport.Children.Remove(lookDirectionRay);
                    lookDirectionRay = null;
                }
                if (lookDirectionRay == null)
                {
                    Point3D[] twopoints = calculateLookVectorPositionStartingFromCurrentPosition();
                    Model3DGroup TheLine = ReturnModel3DGroupPerLine(5, twopoints[0], twopoints[1], System.Windows.Media.Brushes.Lavender);
                    lookDirectionRay = new ModelVisual3D() { Content = TheLine };
                    viewport.Children.Add(lookDirectionRay);
                }
                
            }
        }
        private Point3D[] calculateLookVectorPositionStartingFromCurrentPosition()
        {
            //angleaccumalated_delet = angleaccumalated_delet + XMovement_InDegress;
            double maxdist = Mycamera.FarPlaneDistance;
            double stam = 0;
            if (Math.Abs(angleaccumalated_delet) == 89 || Math.Abs(angleaccumalated_delet) == 90)
            {
                stam = stam + 1;
            }
            stam = stam + 1;
            Vector3D currentUpDirection = Mycamera.UpDirection;
            Vector3D currentLookDirection = Mycamera.LookDirection;
            Point3D currentPosition = Mycamera.Position;
            Vector3D currentLookDirection_untivaector = currentLookDirection;
            currentLookDirection_untivaector.Normalize();
            Point3D p1 = currentPosition + currentLookDirection_untivaector * 200;
            Vector3D distVector;
            Vector3D normVectorZXPlane = new Vector3D(0, 1, 0);
            Vector3D normVectorYXPlane = new Vector3D(0, 0, 1);
            Vector3D normVectorYZPlane = new Vector3D(1, 0, 0);
            Point3D pointLookedAtConstiutingRottaionCenter;
            bool res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYXPlane, out pointLookedAtConstiutingRottaionCenter);
            if (res == false)
            {
                res = VectorProjections.findIntersectionOfVectorOnPalneStatringAtPointP(currentPosition, currentLookDirection, normVectorYZPlane, out pointLookedAtConstiutingRottaionCenter);
                //normVectorZXPlane = normVectorYZPlane;
            }
            Vector3D lookVectorTowardPlaneInQuestion = -(currentPosition - pointLookedAtConstiutingRottaionCenter);
            Vector3D Proj;
            if (findDistanceBetween2Vector(new Vector3D(1, 0, 0), currentUpDirection) < 0.05)
            {
                Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorYZPlane, out  distVector);//
            }
            else
            {
                Proj = VectorProjections.calculateProjectionOfVectorIntoPlane(lookVectorTowardPlaneInQuestion, normVectorZXPlane, out  distVector);//
            }
            Point3D p2 = p1 + currentLookDirection_untivaector * maxdist;
            return new Point3D[] { p1,p2};
        }
        void CompositionTarget_Rendering_V2(object sender, EventArgs e)
        {
            //return;
            //ConsoleManager.Show();
            //Console.WriteLine("hello");
            if (lastRenderTime == ((RenderingEventArgs)e).RenderingTime)
                return;

            lastRenderTime = ((RenderingEventArgs)e).RenderingTime;
            CompositionTarget_Rendering_V2_inaction();
            //displayRayIfNecessary();

        }
        void CompositionTarget_Rendering_V2_backup(object sender, EventArgs e)
        {
            //return;
            //ConsoleManager.Show();
            //Console.WriteLine("hello");
            if (lastRenderTime == ((RenderingEventArgs)e).RenderingTime)
                return;

            lastRenderTime = ((RenderingEventArgs)e).RenderingTime;
            foreach (TubeGeometry TG in ShapesList)
            {







                TG.ExternalAxis_v2.CalculateMatrix();
                //TG.ExternalAxis.ApplyMAtrixToDixtionaryOfPoints(TG.Important3DPoints);

                Matrix3D RotationAroundExternalAxis = TG.ExternalAxis_v2.RotationAroundArbAxix.RotationMatrix3D;

                ((Transform3DGroup)(TG.GeometriesAndMtrices.mGeometry.Transform)).Children[0] = new MatrixTransform3D(RotationAroundExternalAxis);
                if (TG.DisplayPole == true)
                {
                    ((Transform3DGroup)(TG.GeometriesAndMtrices.mGeometry3.Transform)).Children[0] = new MatrixTransform3D(RotationAroundExternalAxis);
                }
                
                if (TG.DisplayCommet == true && 1==0)
                {
                    Matrix3D TranslateMatrixInitial2 = returnMatrixPertransformation(0, TG.TheAxis, new Vector3D(0, 0, 0));//selfrotatingangle was removed
                    ((Transform3DGroup)(TG.GeometriesAndMtrices.mGeometry2.Transform)).Children[0] = new MatrixTransform3D(TranslateMatrixInitial2);
                }
                if (TG.DisplayTrailFlag == true)
                {
                    TG.Trail.AddPointIfApplicable(TG.Important3DPoints_V2["CENTER POINT"].Now);
                }

                //SpinMode = SpinMode + 1;
                //if (SpinMode == 4) SpinMode = 0;
                //break;
            }


        }
        void CompositionTarget_Rendering(object sender, EventArgs e)
        {

            if (lastRenderTime == ((RenderingEventArgs)e).RenderingTime)
                return;

            lastRenderTime = ((RenderingEventArgs)e).RenderingTime;
            foreach (TubeGeometry TG in ShapesList)
            {


                //TG.selfrotatingangle = TG.selfrotatingangle + .8;

                //if (((TG as Ellipsoid) == null))
                //{
                //    TG.groundrotatingangle = TG.groundrotatingangle + .4;
                //}


                //RecalculateTheAxis(TG.groundrotatingangle);





                TG.ExternalAxis.CalculateMatrix();
                TG.ExternalAxis.ApplyMAtrixToDixtionaryOfPoints(TG.Important3DPoints);

                Matrix3D RotationAroundExternalAxis = TG.ExternalAxis.RotationAroundArbAxix.RotationMatrix3D;

                ((Transform3DGroup)(TG.GeometriesAndMtrices.mGeometry.Transform)).Children[0] = new MatrixTransform3D(RotationAroundExternalAxis);

                if (TG.DisplayCommet == true && 1 == 0)
                {
                    Matrix3D TranslateMatrixInitial2 = returnMatrixPertransformation(0, TG.TheAxis, new Vector3D(0, 0, 0));//selfrotatingangle was removed
                    ((Transform3DGroup)(TG.GeometriesAndMtrices.mGeometry2.Transform)).Children[0] = new MatrixTransform3D(TranslateMatrixInitial2);
                }


                //SpinMode = SpinMode + 1;
                //if (SpinMode == 4) SpinMode = 0;
            }


        }
        void CompositionTarget_Rendering_v1(object sender, EventArgs e)
        {
            //TubeGeometry TG;
            //TG = ShapesList[0];

                //testmatrix();
                ////return;
                //// Ensure we only do this once per frame
                //Vector3D TestVector = new Vector3D(4, 3, 0);
                //Vector3D RestlVectoTest;
                //Vector3D RestlVectoTest2;
                //Vector3D RestlVectoTest3;
                if (lastRenderTime == ((RenderingEventArgs)e).RenderingTime)
                    return;

                lastRenderTime = ((RenderingEventArgs)e).RenderingTime;
                foreach (TubeGeometry TG in ShapesList)
                {
                //rotationY.Angle = rotationY.Angle + 1;
                //if (SpinMode <= 2)
                //{
                //    selfrotatingangle = selfrotatingangle + 2;
                //}
                //else
                //{
                //    groundrotatingangle = groundrotatingangle + 2;
                //}

                TG.selfrotatingangle = TG.selfrotatingangle + .8;
                //TG.selfrotatingangle = 0;
                if (((TG as Ellipsoid) == null))
                {
                    TG.groundrotatingangle = TG.groundrotatingangle + .4;
                }
                //TG.selfrotatingangle = 90;

                RecalculateTheAxis(TG.groundrotatingangle);

               
                Matrix3D TranslateMatrixInitial = returnMatrixPertransformation(0, TG.TheAxis, -TG.TranslationVector);//selfrotatingangle was removed
                Matrix3D WheelLikeSpinMatrix = returnMatrixPertransformation(TG.selfrotatingangle, TG.TheAxis, new Vector3D(0, 0, 0));//selfrotatingangle was removed
                Matrix3D GroundRotataionMatrix = returnMatrixPertransformation(TG.groundrotatingangle, TG.GroundAxis, new Vector3D(0, 0, 0));//groundrotatingangle  was removed
                Matrix3D TranslateBackMatrix = returnMatrixPertransformation(0, TG.GroundAxis, TG.TranslationVector);//selfrotatingangle was removed


                Matrix3D TranslateBackMatrixTest = returnMatrixPertransformation(0, TG.GroundAxis, -TG.TranslationVector);
                Matrix3D WheelLikeSpinMatrixTest = returnMatrixPertransformation(TG.selfrotatingangle, TG.TheAxis, new Vector3D(0, 0, 0));//selfrotatingangle was removed
                //RestlVectoTest = Vector3D.Multiply(TestVector, TranslateBackMatrixTest);
                //RestlVectoTest2 = Vector3D.Multiply(RestlVectoTest, WheelLikeSpinMatrixTest);
                //RestlVectoTest3 = Vector3D.Multiply(RestlVectoTest2, GroundRotataionMatrix); 
                //int t = 0;
                //RestlVectoTest2 = Vector3D.Multiply(RestlVectoTest, GroundRotataionMatrix); ;
                //Matrix3D retmtarix3D = Matrix3D.Multiply(matrix1, matrix2);
                //matrix2.Append(matrix1);
                //GroundRotataionMatrix.Append(TranslateBackMatrix);
                TranslateMatrixInitial.Append(WheelLikeSpinMatrix);
                GroundRotataionMatrix.Append(TranslateBackMatrix);
                TranslateMatrixInitial.Append(GroundRotataionMatrix);
                //WheelLikeSpinMatrix.Append(TranslateBackMatrix);
                //RestlVectoTest3 = Vector3D.Multiply(TestVector, WheelLikeSpinMatrix); ;
                //mGeometry.Transform.Value = matrix2;
                //myMatrixTransform3D = matrix2;
                ((Transform3DGroup)(TG.GeometriesAndMtrices.mGeometry.Transform)).Children[0] = new MatrixTransform3D(TranslateMatrixInitial);
                if (TG.DisplayCommet == true)
                {
                    Matrix3D TranslateMatrixInitial2 = returnMatrixPertransformation(0, TG.TheAxis, new Vector3D(0, 0, 0));//selfrotatingangle was removed
                    ((Transform3DGroup)(TG.GeometriesAndMtrices.mGeometry2.Transform)).Children[0] = new MatrixTransform3D(TranslateMatrixInitial2);
                }


                SpinMode = SpinMode + 1;
                if (SpinMode == 4) SpinMode = 0;
            }


        }

       //SHapeParameters
        public SHapeParameters ReturnShape5()
        {
            
            Bitmap ResourceBitmap = (Bitmap)Wpf3DTube_NonStatic.Properties.Resources.IMG_0804x;

            BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
                                              ResourceBitmap.GetHbitmap(),
                                              IntPtr.Zero,
                                              System.Windows.Int32Rect.Empty,
                                              BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            //ImageBrush ib2 = new ImageBrush(bs); 
            //ImageBrush ib2 = ReturnImageBrusgBasedOnFileName();

            //object ibx = TubeGeometry.CreateCheckersDrawingBrush();
            object ibx = TubeGeometry.CreateTransparemtFramedBrush();
           
            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = new Vector3D(0, 1, -.2);
            TranslationVector = 1.25*new Vector3D(10, -2, 5);

            double TubeThickness = 3;
            double InternalDiameter = 8;
            int angleStep_horizotal = 4;
            int angleStep_vertical = 10;
            int NumberOfDefreesToComplete_horizontal = 360;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 0;
            double intialGroundrotation = 0;
            bool isitTile = true;
            bool displaycommet = false;
            Type MyType = typeof(TubeGeometry);
            double AngleMovementAroundExternalAxis = 1;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = 100;
            double InitialOrbitalAngle = 0;
            double Tilt = 25;
            bool DisplayAxixPole = false;
            bool DisplayTrail = true;
            PlanetPropeties PlanetRelated = new PlanetPropeties();
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ibx,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated
                                         ,0
                                         ,0
                                         ,0,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape2()
        {
            Bitmap ResourceBitmap = (Bitmap)Wpf3DTube_NonStatic.Properties.Resources.IMG_0804x;

            BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
                                              ResourceBitmap.GetHbitmap(),
                                              IntPtr.Zero,
                                              System.Windows.Int32Rect.Empty,
                                              BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            //ImageBrush ib2 = new ImageBrush(bs); 
            ImageBrush ib2 = ReturnImageBrusgBasedOnFileName();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = new Vector3D(0, 1, 0);
            TranslationVector =( new Vector3D(-5, -8, 5))*1.5;

            double TubeThickness = 2;
            double InternalDiameter = 5;
            int angleStep_horizotal = 2;
            int angleStep_vertical = 2;
            int NumberOfDefreesToComplete_horizontal = 360;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 90;
            double intialGroundrotation = 90;
            bool isitTile = false;
            bool displaycommet = false;
            Type MyType = typeof(TubeGeometry);
            double AngleMovementAroundExternalAxis = 1;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = 0;
            double InitialOrbitalAngle = 0;
            double Tilt = 0;
            bool DisplayAxixPole = false;
            bool DisplayTrail = false;
            PlanetPropeties PlanetRelated = new PlanetPropeties();
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ib2,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated
                                         ,0
                                         ,0
                                         ,0,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape3()
        {
            Bitmap ResourceBitmap = (Bitmap)Wpf3DTube_NonStatic.Properties.Resources.IMG_0804x;

            BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
                                              ResourceBitmap.GetHbitmap(),
                                              IntPtr.Zero,
                                              System.Windows.Int32Rect.Empty,
                                              BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            //ImageBrush ib2 = new ImageBrush(bs); 
            ImageBrush ib2 = ReturnImageBrusgBasedOnFileName();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = new Vector3D(0, 1, -.3);
             
            TranslationVector = 1*(new Vector3D(-15, 0, -15));

            double TubeThickness = 2;
            double InternalDiameter = 6;
            int angleStep_horizotal = 2;
            int angleStep_vertical = 2;
            int NumberOfDefreesToComplete_horizontal = 360;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 45;
            double intialGroundrotation = 45;
            bool isitTile = false;
            bool displaycommet = false;
            Type MyType = typeof(TubeGeometry);
            double AngleMovementAroundExternalAxis = 1;
            double AngleMovementAroundInternalAxis = 1;
            double  DistanceFromRespectiveSun = 75;
            double InitialOrbitalAngle = 0;
            double Tilt = 0;
            bool DisplayAxixPole = false;
            bool DisplayTrail = true;
            PlanetPropeties PlanetRelated = new PlanetPropeties();
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ib2,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated
                                         ,0
                                         ,0
                                         ,0,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        
                                        );
        }


        public SHapeParameters ReturnShape4()
        {
            //Bitmap ResourceBitmap = (Bitmap)Wpf3DTube_NonStatic.Properties.Resources.IMG_0804x;

            //BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
            //                                  ResourceBitmap.GetHbitmap(),
            //                                  IntPtr.Zero,
            //                                  System.Windows.Int32Rect.Empty,
            //                                  BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            ////ImageBrush ib2 = new ImageBrush(bs); 
            //ImageBrush ib2 = ReturnImageBrusgBasedOnFileName();
            object ibx = TubeGeometry.CreateTransparemtFramedBrush();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = new Vector3D(0,1,0);
            //TranslationVector = 1 * (new Vector3D(0,0,0));
            TranslationVector = 1 * (new Vector3D(-10, 5, -5));


            int angleStep_horizotal = 10;
            int angleStep_vertical = 3;
            int NumberOfDefreesToComplete_horizontal = 180;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 45;
            double intialGroundrotation = 45;
            bool isitTile = true;
            bool displaycommet = false;
            Type MyType = typeof(Ellipsoid);

            //double EllipsisA = 15;
            //double EllipsisB = 6;
            //double EllipsisC = 10;
            double EllipsisA = 10;
            double EllipsisB = 10;
            double EllipsisC = 10;

            double TubeThickness = 1;
            double InternalDiameter = 1 ;
            double AngleMovementAroundExternalAxis = 1;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = 100;
            double InitialOrbitalAngle = 0;
            double Tilt = 0;
            bool DisplayAxixPole = false;
            bool DisplayTrail = false;
            PlanetPropeties PlanetRelated = new PlanetPropeties();
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ibx,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated,
                                         EllipsisA,
                                         EllipsisB,
                                         EllipsisC,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape1_v2(Point3D center)
        {
            //Bitmap ResourceBitmap = (Bitmap)Wpf3DTube_NonStatic.Properties.Resources.IMG_0804x;

            //BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
            //                                  ResourceBitmap.GetHbitmap(),
            //                                  IntPtr.Zero,
            //                                  System.Windows.Int32Rect.Empty,
            //                                  BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            ////ImageBrush ib2 = new ImageBrush(bs); 
            //ImageBrush ib2 = ReturnImageBrusgBasedOnFileName();
            object ibx = TubeGeometry.CreateTransparemtFramedBrush();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = new Vector3D(-.3, 1, 0);
            //TranslationVector = 1 * (new Vector3D(0,0,0));
            TranslationVector = (Vector3D)center;// new Vector3D(-10, 5, -5);


            int angleStep_horizotal = 10;
            int angleStep_vertical = 3;
            int NumberOfDefreesToComplete_horizontal = 180;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 45;
            double intialGroundrotation = 45;
            bool isitTile = true;
            bool displaycommet = false;
            Type MyType = typeof(Ellipsoid);

            //double EllipsisA = 15;
            //double EllipsisB = 6;
            //double EllipsisC = 10;
            double EllipsisA = 10;
            double EllipsisB = 10;
            double EllipsisC = 10;

            double TubeThickness = 1;
            double InternalDiameter = 1;
            double AngleMovementAroundExternalAxis = 1;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = 70;
            double InitialOrbitalAngle = 0;
            double Tilt = 0;
            bool DisplayAxixPole = false;
            bool DisplayTrail = false;
            PlanetPropeties PlanetRelated = new PlanetPropeties();
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ibx,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated,
                                         EllipsisA,
                                         EllipsisB,
                                         EllipsisC,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape1(Point3D center)
        {
            //Bitmap ResourceBitmap = (Bitmap)Wpf3DTube_NonStatic.Properties.Resources.IMG_0804x;

            //BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
            //                                  ResourceBitmap.GetHbitmap(),
            //                                  IntPtr.Zero,
            //                                  System.Windows.Int32Rect.Empty,
            //                                  BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            ////ImageBrush ib2 = new ImageBrush(bs); 
            //ImageBrush ib2 = ReturnImageBrusgBasedOnFileName();
            object ibx = TubeGeometry.CreateTransparemtFramedBrush();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = new Vector3D(-.3, 1, 0);
            //TranslationVector = 1 * (new Vector3D(0,0,0));
            TranslationVector = (Vector3D)center;// new Vector3D(-10, 5, -5);


            int angleStep_horizotal = 10;
            int angleStep_vertical = 3;
            int NumberOfDefreesToComplete_horizontal = 180;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 45;
            double intialGroundrotation = 45;
            bool isitTile = true;
            bool displaycommet = false;
            Type MyType = typeof(Ellipsoid);

            //double EllipsisA = 15;
            //double EllipsisB = 6;
            //double EllipsisC = 10;
            double EllipsisA = 10;
            double EllipsisB = 10;
            double EllipsisC = 10;

            double TubeThickness = 1;
            double InternalDiameter = 1;
            double AngleMovementAroundExternalAxis = 1;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = 70;
            double InitialOrbitalAngle = 0;
            double Tilt = 0;
            bool DisplayAxixPole = false;
            bool DisplayTrail = false;
            PlanetPropeties PlanetRelated = new PlanetPropeties();
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ibx,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated,
                                         EllipsisA,
                                         EllipsisB,
                                         EllipsisC,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape6()
        {
            Bitmap ResourceBitmap = (Bitmap)Wpf3DTube_NonStatic.Properties.Resources.IMG_0804x;

            BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
                                              ResourceBitmap.GetHbitmap(),
                                              IntPtr.Zero,
                                              System.Windows.Int32Rect.Empty,
                                              BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            //ImageBrush ib2 = new ImageBrush(bs); 
            //ImageBrush ib2 = ReturnImageBrusgBasedOnFileName();

            //object ibx = TubeGeometry.CreateCheckersDrawingBrush();
            object ibx = TubeGeometry.CreateTransparemtFramedBrush();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = new Vector3D(-.3, 1, 0);
            TranslationVector = 1.25 * new Vector3D(10, -2, 5);

            double TubeThickness = .9;
            double InternalDiameter = 25;
            int angleStep_horizotal = 4;
            int angleStep_vertical = 10;
            int NumberOfDefreesToComplete_horizontal = 360;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 0;
            double intialGroundrotation = 0;
            bool isitTile = true;
            bool displaycommet = false;
            Type MyType = typeof(TubeGeometry);
            double AngleMovementAroundExternalAxis = 1;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = 70;
            double InitialOrbitalAngle = 0;
            double Tilt = 0;
            bool DisplayAxixPole = false;
            bool DisplayTrail = true;
            PlanetPropeties PlanetRelated = new PlanetPropeties();
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ibx,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated
                                         , 0
                                         , 0
                                         , 0,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape6_v2()
        {
            PlanetPropeties earth = Planets["earth"];
            object ibx = TubeGeometry.CreateTransparemtFramedBrush();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = earth.mormalToPlanetEquator;// new Vector3D(0, 1, 0);
            AxisVector = new Vector3D(0, 1, 0);
            //TranslationVector = 1 * (new Vector3D(0,0,0));
            TranslationVector = (Vector3D)earth.initilaLocation;// 1 * (new Vector3D(-10, 5, -5));
            //TranslationVector = new Vector3D();

            int angleStep_horizotal = 10;
            int angleStep_vertical = 3;
            int NumberOfDefreesToComplete_horizontal = 180;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 0;
            double intialGroundrotation = 45;
            bool isitTile = true;
            bool displaycommet = false;
            Type MyType = typeof(Ellipsoid);

            double EllipsisA = earth.Dimensions.EquaterialSemiMajorLength / AstronomicalUnint;
            double EllipsisB =1* earth.Dimensions.EquaterialSemiMinorLength / AstronomicalUnint;
            double EllipsisC = 1*earth.Dimensions.PolarRadius / AstronomicalUnint;

            double TubeThickness = 1;
            double InternalDiameter = 1;
            double AngleMovementAroundExternalAxis = 0.5;// 1.0 / 20.0;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = earth.OrbitSemiMinorAxis / AstronomicalUnint;
            double InitialOrbitalAngle = 0;
            //double Tilt = 0;// earth.Tilt;
            double Tilt =  earth.Tilt;
            bool DisplayAxixPole = true;
            bool DisplayTrail = true;
            PlanetPropeties PlanetRelated = earth;
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ibx,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated,
                                         EllipsisA,
                                         EllipsisB,
                                         EllipsisC,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape6_v2_bkup()
        {
            PlanetPropeties earth = Planets["earth"];
            object ibx = TubeGeometry.CreateTransparemtFramedBrush();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = earth.mormalToPlanetEquator;// new Vector3D(0, 1, 0);
            //TranslationVector = 1 * (new Vector3D(0,0,0));
            TranslationVector = (Vector3D)earth.initilaLocation;// 1 * (new Vector3D(-10, 5, -5));


            int angleStep_horizotal = 10;
            int angleStep_vertical = 3;
            int NumberOfDefreesToComplete_horizontal = 180;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 0;
            double intialGroundrotation = 45;
            bool isitTile = true;
            bool displaycommet = false;
            Type MyType = typeof(Ellipsoid);

            double EllipsisA = earth.Dimensions.EquaterialSemiMajorLength / AstronomicalUnint;
            double EllipsisB = earth.Dimensions.EquaterialSemiMinorLength/AstronomicalUnint;
            double EllipsisC = earth.Dimensions.PolarRadius / AstronomicalUnint;

            double TubeThickness = 1;
            double InternalDiameter = 1;
            double AngleMovementAroundExternalAxis = 1.0/5.0;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = earth.OrbitSemiMinorAxis / AstronomicalUnint;
            double InitialOrbitalAngle = 0;
            double Tilt = earth.Tilt;
            bool DisplayAxixPole = true;
            bool DisplayTrail = true;
            PlanetPropeties PlanetRelated = earth;
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ibx,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated,
                                         EllipsisA,
                                         EllipsisB,
                                         EllipsisC,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape61_v2()
        {
            PlanetPropeties earth = Planets["earth2"];
            object ibx = TubeGeometry.CreateTransparemtFramedBrush();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = earth.mormalToPlanetEquator;// new Vector3D(0, 1, 0);
            //TranslationVector = 1 * (new Vector3D(0,0,0));
            TranslationVector = (Vector3D)earth.initilaLocation;// 1 * (new Vector3D(-10, 5, -5));


            int angleStep_horizotal = 10;
            int angleStep_vertical = 3;
            int NumberOfDefreesToComplete_horizontal = 180;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 0;
            double intialGroundrotation = 45;
            bool isitTile = true;
            bool displaycommet = false;
            Type MyType = typeof(Ellipsoid);

            double EllipsisA = earth.Dimensions.EquaterialSemiMajorLength / AstronomicalUnint;
            double EllipsisB = earth.Dimensions.EquaterialSemiMinorLength / AstronomicalUnint;
            double EllipsisC = earth.Dimensions.PolarRadius / AstronomicalUnint;

            double TubeThickness = 1;
            double InternalDiameter = 1;
            double AngleMovementAroundExternalAxis = 1.0 / 5.0;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = earth.OrbitSemiMinorAxis / AstronomicalUnint;
            double InitialOrbitalAngle = 0;
            double Tilt = earth.Tilt;
            bool DisplayAxixPole = true;
            bool DisplayTrail = true;
            PlanetPropeties PlanetRelated = earth;
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ibx,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated,
                                         EllipsisA,
                                         EllipsisB,
                                         EllipsisC,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape7()
        {
            Bitmap ResourceBitmap = (Bitmap)Wpf3DTube_NonStatic.Properties.Resources.IMG_0804x;

            BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
                                              ResourceBitmap.GetHbitmap(),
                                              IntPtr.Zero,
                                              System.Windows.Int32Rect.Empty,
                                              BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            //ImageBrush ib2 = new ImageBrush(bs); 
            ImageBrush ib2 = ReturnImageBrusgBasedOnFileName();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = new Vector3D(0, 1, 0);
            TranslationVector = (new Vector3D(-5, -8, 5)) * 1.5;

            double TubeThickness = 10;
            double InternalDiameter = 12;
            int angleStep_horizotal = 2;
            int angleStep_vertical = 2;
            int NumberOfDefreesToComplete_horizontal = 360;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 90;
            double intialGroundrotation = 90;
            bool isitTile = false;
            bool displaycommet = false;
            Type MyType = typeof(TubeGeometry);
            double AngleMovementAroundExternalAxis = 1;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = 80;
            double InitialOrbitalAngle = 0;
            double Tilt = 0;
            bool DisplayAxixPole = false;
            bool DisplayTrail = false;
            PlanetPropeties PlanetRelated = new PlanetPropeties();
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ib2,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated
                                         , 0
                                         , 0
                                         , 0,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape8()
        {
            Bitmap ResourceBitmap = (Bitmap)Wpf3DTube_NonStatic.Properties.Resources.IMG_0804x;

            BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
                                              ResourceBitmap.GetHbitmap(),
                                              IntPtr.Zero,
                                              System.Windows.Int32Rect.Empty,
                                              BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            //ImageBrush ib2 = new ImageBrush(bs); 
            ImageBrush ib2 = ReturnImageBrusgBasedOnFileName();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = new Vector3D(0, 1, 0);
            TranslationVector = (new Vector3D(-5, -8, 5)) * 1.5;

            double TubeThickness = 5;
            double InternalDiameter = 3;
            int angleStep_horizotal = 2;
            int angleStep_vertical = 2;
            int NumberOfDefreesToComplete_horizontal = 360;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 90;
            double intialGroundrotation = 90;
            bool isitTile = false;
            bool displaycommet = false;
            Type MyType = typeof(TubeGeometry);
            double AngleMovementAroundExternalAxis = 1;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = 80;
            double InitialOrbitalAngle = 0;
            double Tilt = 15;
            bool DisplayAxixPole = false;
            bool DisplayTrail = false;
            PlanetPropeties PlanetRelated = new PlanetPropeties();
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ib2,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated
                                         , 0
                                         , 0
                                         , 0,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape9_v2()
        {

            PlanetPropeties sun = Planets["sun"];
            object ibx = TubeGeometry.CreateTransparemtFramedBrush();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = sun.mormalToPlanetEquator;// new Vector3D(0, 1, 0);
            //TranslationVector = 1 * (new Vector3D(0,0,0));
            TranslationVector = (Vector3D)sun.initilaLocation;// 1 * (new Vector3D(-10, 5, -5));


            int angleStep_horizotal = 10;
            int angleStep_vertical = 3;
            int NumberOfDefreesToComplete_horizontal = 180;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 0;
            double intialGroundrotation = 45;
            bool isitTile = true;
            bool displaycommet = false;
            Type MyType = typeof(Ellipsoid);

            //double EllipsisA = 15;
            //double EllipsisB = 6;
            //double EllipsisC = 10;
            double EllipsisA = sun.Dimensions.EquaterialSemiMajorLength / AstronomicalUnint;
            double EllipsisB = sun.Dimensions.EquaterialSemiMinorLength / AstronomicalUnint;
            double EllipsisC = sun.Dimensions.PolarRadius / AstronomicalUnint;

            double TubeThickness = 1;
            double InternalDiameter = 1;
            double AngleMovementAroundExternalAxis = 1;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = 0;
            double InitialOrbitalAngle = 0;
            double Tilt = 0;
            bool DisplayAxixPole = false;
            bool DisplayTrail = false;
            PlanetPropeties PlanetRelated = sun;
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ibx,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated,
                                         EllipsisA,
                                         EllipsisB,
                                         EllipsisC,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }
        public SHapeParameters ReturnShape9()
        {
            //Bitmap ResourceBitmap = (Bitmap)Wpf3DTube_NonStatic.Properties.Resources.IMG_0804x;

            //BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
            //                                  ResourceBitmap.GetHbitmap(),
            //                                  IntPtr.Zero,
            //                                  System.Windows.Int32Rect.Empty,
            //                                  BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            ////ImageBrush ib2 = new ImageBrush(bs); 
            //ImageBrush ib2 = ReturnImageBrusgBasedOnFileName();
            object ibx = TubeGeometry.CreateTransparemtFramedBrush();

            Vector3D AxisVector;
            Vector3D TranslationVector;
            AxisVector = new Vector3D(0, 1, 0);
            //TranslationVector = 1 * (new Vector3D(0,0,0));
            TranslationVector = 1 * (new Vector3D(-10, 5, -5));


            int angleStep_horizotal = 10;
            int angleStep_vertical = 3;
            int NumberOfDefreesToComplete_horizontal = 180;
            int NumberOfDefreesToComplete_vertical = 360;
            double intialSelfRotation = 45;
            double intialGroundrotation = 45;
            bool isitTile = true;
            bool displaycommet = false;
            Type MyType = typeof(Ellipsoid);

            //double EllipsisA = 15;
            //double EllipsisB = 6;
            //double EllipsisC = 10;
            double EllipsisA = 40;
            double EllipsisB = 40;
            double EllipsisC = 40;

            double TubeThickness = 1;
            double InternalDiameter = 1;
            double AngleMovementAroundExternalAxis = 1;
            double AngleMovementAroundInternalAxis = 1;
            double DistanceFromRespectiveSun = 0;
            double InitialOrbitalAngle = 0;
            double Tilt = 0;
            bool DisplayAxixPole = false;
            bool DisplayTrail = false;
            PlanetPropeties PlanetRelated = new PlanetPropeties();
            return new SHapeParameters(
                                         AxisVector,
                                         TranslationVector,
                                         TubeThickness,
                                         InternalDiameter,
                                         angleStep_horizotal,
                                         angleStep_vertical,
                                         NumberOfDefreesToComplete_horizontal,
                                         NumberOfDefreesToComplete_vertical,
                                         ibx,
                                         intialSelfRotation,
                                         intialGroundrotation,
                                         isitTile,
                                         displaycommet,
                                         MyType,
                                         PlanetRelated,
                                         EllipsisA,
                                         EllipsisB,
                                         EllipsisC,
                                         AngleMovementAroundExternalAxis,
                                         AngleMovementAroundInternalAxis,
                                         DistanceFromRespectiveSun,
                                         InitialOrbitalAngle,
                                         Tilt,
                                         DisplayAxixPole,
                                         DisplayTrail

                                        );
        }

        public List<SHapeParameters> ReturnSHapesList_v2(Point3D center)
        {
            List<SHapeParameters> ListShapesParameters = new List<SHapeParameters>();
            //ListShapesParameters.Add(ReturnShape2());
            ListShapesParameters.Add(ReturnShape9_v2()); //sun





            ListShapesParameters.Add(ReturnShape6_v2_bkup()); // orbitting
            //ListShapesParameters.Add(ReturnShape61_v2()); // orbitting
            
            //ListShapesParameters.Add(ReturnShape1_v2(center));






            ////ListShapesParameters.Add(ReturnShape7());
            //ListShapesParameters.Add(ReturnShape8());

            //ListShapesParameters.Add(ReturnShape3());
            //ListShapesParameters.Add(ReturnShape4());
            //ListShapesParameters.Add(ReturnShape5());



            return ListShapesParameters;
        }
        public List<SHapeParameters> ReturnSHapesList(Point3D center)
        {
            List<SHapeParameters> ListShapesParameters = new List<SHapeParameters>();
            //ListShapesParameters.Add(ReturnShape2());
            ListShapesParameters.Add(ReturnShape9());





            ListShapesParameters.Add(ReturnShape6());
            ListShapesParameters.Add(ReturnShape1(center));






            ////ListShapesParameters.Add(ReturnShape7());
            //ListShapesParameters.Add(ReturnShape8());

            ListShapesParameters.Add(ReturnShape3());
            ListShapesParameters.Add(ReturnShape4());
            ListShapesParameters.Add(ReturnShape5());

            
            
            return ListShapesParameters;
        }

        public void Build3DTube_v2(List<SHapeParameters> Shapes, Point3D center,Vector3D axisToAll)
        {

            
            AxisAndCenter AxisAndCenterExternal;
            AxisAndCenter AxisAndCenterExternal_v2;
            double SinRotorOverCosRotor = 1;// in the future this will be a parameter per each shape
            double ExternalAngleMovement;
            double InternalAngleMovement;
            double ratioBetwenAxis;
            double prTiltRelatedRadius = 0;
            int counter = 0;
            foreach (SHapeParameters shape1 in Shapes)
            {

                center = (Point3D)(((Vector3D)shape1.prPlanetRelated.imaginaryCenter) / AstronomicalUnint);
                double ellipsisA; double ellipsisB; double ellipsisC;
                ellipsisA = shape1.prEllipsisA;
                ellipsisB = shape1.prEllipsisB;
                ellipsisC = shape1.prEllipsisC;
                ExternalAngleMovement = shape1.prAngleMovementAroundExternalAxis;
                InternalAngleMovement = shape1.prAngleMovementAroundInternalAxis;
                AxisAndCenterExternal = new AxisAndCenter(shape1.prTilt, ExternalAngleMovement, axisToAll,InternalAngleMovement, center, SinRotorOverCosRotor,prTiltRelatedRadius);
                counter++;
                Type EllipisoidType = typeof(Ellipsoid);
                TubeGeometry TG;
                if (shape1.prShapeType == EllipisoidType)
                {
                    prTiltRelatedRadius = ellipsisA / 2;
                    //prTiltRelatedRadius = ellipsisA ;
                }
                else
                {
                    prTiltRelatedRadius = shape1.prTubeThickness/2;
                }
                if (counter==1) 
                {
                    ratioBetwenAxis = SinRotorOverCosRotor*1;
                } else 
                {
                    ratioBetwenAxis = SinRotorOverCosRotor;
                }
                if (shape1.prPlanetRelated.OrbitSemiMajorAxis == 0)
                {
                    ratioBetwenAxis = 1;
                }
                else
                {
                    ratioBetwenAxis =  shape1.prPlanetRelated.OrbitSemiMinorAxis / shape1.prPlanetRelated.OrbitSemiMajorAxis;
                }
                    AxisAndCenterExternal_v2 = new AxisAndCenter_v2(shape1.prTilt, ExternalAngleMovement, shape1.prAxisVector, InternalAngleMovement, center, ratioBetwenAxis, prTiltRelatedRadius);
                Model3DGroup ModelDgroup = ReturnModel3DGroup();
                Build3DTube_inaction_v2(ModelDgroup, shape1.prTranslationVector, shape1.prImageBrush, shape1.prAxisVector, shape1.prTubeThickness, shape1.prInternalDiameter, shape1.prangleStep_horizotal, shape1.prangleStep_vertical, shape1.prNumberOfDefreesToComplete_horizontal, shape1.prNumberOfDefreesToComplete_vertical, shape1.prinitialselfRotation, shape1.prinitialGroundAngle, shape1.prIsItTile, shape1.prDisplayCommet, shape1.prShapeType, ellipsisA, ellipsisB, ellipsisC, AxisAndCenterExternal, AxisAndCenterExternal_v2, ExternalAngleMovement, InternalAngleMovement, shape1.prDistanceFromRespectiveSun, shape1.prInitialOrbitalAngle, shape1.prTilt, shape1.prDisplayAxixPole, shape1.prDisplayTrail, shape1.prPlanetRelated);
            }
        }
        public void Build3DTube(List<SHapeParameters> Shapes, Point3D center, Vector3D axisToAll)
        {

            AxisAndCenter AxisAndCenterExternal;
            AxisAndCenter AxisAndCenterExternal_v2;
            double SinRotorOverCosRotor = 1;// in the future this will be a parameter per each shape
            double ExternalAngleMovement;
            double InternalAngleMovement;
            double ratioBetwenAxis;
            double prTiltRelatedRadius = 0;
            int counter = 0;
            foreach (SHapeParameters shape1 in Shapes)
            {

                double ellipsisA; double ellipsisB; double ellipsisC;
                ellipsisA = shape1.prEllipsisA;
                ellipsisB = shape1.prEllipsisB;
                ellipsisC = shape1.prEllipsisC;
                ExternalAngleMovement = shape1.prAngleMovementAroundExternalAxis;
                InternalAngleMovement = shape1.prAngleMovementAroundInternalAxis;
                AxisAndCenterExternal = new AxisAndCenter(shape1.prTilt, ExternalAngleMovement, axisToAll, InternalAngleMovement, center, SinRotorOverCosRotor, prTiltRelatedRadius);
                counter++;
                Type EllipisoidType = typeof(Ellipsoid);
                TubeGeometry TG;
                if (shape1.prShapeType == EllipisoidType)
                {
                    prTiltRelatedRadius = ellipsisA / 2;
                }
                else
                {
                    prTiltRelatedRadius = shape1.prTubeThickness / 2;
                }
                if (counter == 1)
                {
                    ratioBetwenAxis = SinRotorOverCosRotor * 1;
                }
                else
                {
                    ratioBetwenAxis = SinRotorOverCosRotor;
                }
                AxisAndCenterExternal_v2 = new AxisAndCenter_v2(shape1.prTilt, ExternalAngleMovement, shape1.prAxisVector, InternalAngleMovement, center, ratioBetwenAxis, prTiltRelatedRadius);
                Model3DGroup ModelDgroup = ReturnModel3DGroup();
                Build3DTube_inaction(ModelDgroup, shape1.prTranslationVector, shape1.prImageBrush, shape1.prAxisVector, shape1.prTubeThickness, shape1.prInternalDiameter, shape1.prangleStep_horizotal, shape1.prangleStep_vertical, shape1.prNumberOfDefreesToComplete_horizontal, shape1.prNumberOfDefreesToComplete_vertical, shape1.prinitialselfRotation, shape1.prinitialGroundAngle, shape1.prIsItTile, shape1.prDisplayCommet, shape1.prShapeType, ellipsisA, ellipsisB, ellipsisC, AxisAndCenterExternal, AxisAndCenterExternal_v2, ExternalAngleMovement, InternalAngleMovement, shape1.prDistanceFromRespectiveSun, shape1.prInitialOrbitalAngle, shape1.prTilt, shape1.prDisplayAxixPole, shape1.prDisplayTrail, shape1.prPlanetRelated);
            }
        }
        private ImageBrush ReturnImageBrusgBasedOnFileName()
        {
            Bitmap ResourceBitmap = (Bitmap)Wpf3DTube_NonStatic.Properties.Resources.IMG_0804x;

            BitmapSource bs = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
                                              ResourceBitmap.GetHbitmap(),
                                              IntPtr.Zero,
                                              System.Windows.Int32Rect.Empty,
                                              BitmapSizeOptions.FromWidthAndHeight(ResourceBitmap.Width, ResourceBitmap.Height));
            ImageBrush ib2 = new ImageBrush(bs);
            return ib2;
        }
        public void Build3DTube_inaction_v2(Model3DGroup ModelDgroup, Vector3D TranslationVector, object ib, Vector3D AxisVector, double TubeThickness, double InternalDiameter, int angleStep_horizotal, int angleStep_vertical, int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, double intialselfRotation, double InitialgroundRotation, bool IsItTile, bool displaycommet, Type shapeType, double ellipsisA, double ellipsisB, double ellipsisC, AxisAndCenter AxisAndCenterExternal, AxisAndCenter AxisAndCenterExternal_v2, double ExternalAngleMovement, double InternalAngleMovement, double DistanceFromRespectiveSun, double InitialOrbitalAngle, double Tilt, bool DisplayAxixPole , bool DisplayTrail,PlanetPropeties PP)
        {


            DrawingBrush ibx = new DrawingBrush();
            Type EllipisoidType=typeof(Ellipsoid);
            TubeGeometry TG;
            if (shapeType == EllipisoidType) 
            {
                TG = new Ellipsoid(ModelDgroup,TranslationVector, AxisVector, TubeThickness, InternalDiameter, angleStep_horizotal, angleStep_vertical, NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, intialselfRotation, InitialgroundRotation, IsItTile, displaycommet, AxisAndCenterExternal, AxisAndCenterExternal_v2,DisplayAxixPole,DisplayTrail, ellipsisA, ellipsisB, ellipsisC, ExternalAngleMovement, InternalAngleMovement, DistanceFromRespectiveSun, InitialOrbitalAngle, Tilt,PP);
            }
            else
            {
                TG = new TubeGeometry(ModelDgroup,TranslationVector, AxisVector, TubeThickness, InternalDiameter, angleStep_horizotal, angleStep_vertical, NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, intialselfRotation, InitialgroundRotation, IsItTile, displaycommet, AxisAndCenterExternal, AxisAndCenterExternal_v2,DisplayAxixPole,DisplayTrail, ExternalAngleMovement, InternalAngleMovement, DistanceFromRespectiveSun, InitialOrbitalAngle,PP, Tilt);

            } 
            ShapesList.Add(TG);
            MeshGeometry3D mesh = TG.Mesh;
            MaterialGroup MatGR = new MaterialGroup();
            MaterialGroup CircelMaterialGroup = new MaterialGroup();
            MaterialGroup AxisPoleMaterialGroup = new MaterialGroup();
            //MaterialGroup BAckMatGR = new MaterialGroup();
            if (ib is ImageBrush)
            {
                MatGR.Children.Add(new DiffuseMaterial(ib as ImageBrush));
            }
            else
            {
                MatGR.Children.Add(new DiffuseMaterial(ib as DrawingBrush));
            }
            mGeometry = new GeometryModel3D(mesh, MatGR);
            mGeometry.BackMaterial = MatGR;
            Transform3DGroup MyTransform3DGroup = new Transform3DGroup();
            MatrixTransform3D myMatrixTransform3D = new MatrixTransform3D();


            Transform3DGroup MyTransform3DGroup2 = new Transform3DGroup();
            MatrixTransform3D myMatrixTransform3D2 = new MatrixTransform3D();

            Transform3DGroup MyTransform3DGroup3 = new Transform3DGroup();
            MatrixTransform3D myMatrixTransform3D3 = new MatrixTransform3D();

            MyTransform3DGroup.Children.Add(myMatrixTransform3D);
            MyTransform3DGroup2.Children.Add(myMatrixTransform3D2);
            MyTransform3DGroup3.Children.Add(myMatrixTransform3D3);
            mGeometry.Transform = MyTransform3DGroup;
            AxisPoleMaterialGroup.Children.Add(new DiffuseMaterial(System.Windows.Media.Brushes.Green));
            CircelMaterialGroup.Children.Add(new DiffuseMaterial(System.Windows.Media.Brushes.Green));

            if (displaycommet == true)
            {
                mGeometry2 = new GeometryModel3D(TG.CircleMesh, CircelMaterialGroup);
                mGeometry2.BackMaterial = CircelMaterialGroup;
                mGeometry2.Transform = MyTransform3DGroup2;
                TG.PopulateGeometries(mGeometry, mGeometry2);
                group2.Children.Add(TG.GeometriesAndMtrices.mGeometry2);
            } else
            {
                mGeometry2 = new GeometryModel3D();
                TG.PopulateGeometries(mGeometry, mGeometry2);
            }

            if (DisplayAxixPole == true)
            {
                mGeometry3 = new GeometryModel3D(TG.AxisPoleMesh, AxisPoleMaterialGroup);
                mGeometry3.BackMaterial = AxisPoleMaterialGroup;
                mGeometry3.Transform = MyTransform3DGroup3;
                TG.GeometriesAndMtrices.mGeometry3 = mGeometry3;
                group3.Children.Add(TG.GeometriesAndMtrices.mGeometry3);
            }

           
            



            
            group.Children.Add(TG.GeometriesAndMtrices.mGeometry);


            object ffff = model;
            object zzzz = ffff;
        }
        public void Build3DTube_inaction(Model3DGroup ModelDgroup, Vector3D TranslationVector, object ib, Vector3D AxisVector, double TubeThickness, double InternalDiameter, int angleStep_horizotal, int angleStep_vertical, int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, double intialselfRotation, double InitialgroundRotation, bool IsItTile, bool displaycommet, Type shapeType, double ellipsisA, double ellipsisB, double ellipsisC, AxisAndCenter AxisAndCenterExternal, AxisAndCenter AxisAndCenterExternal_v2, double ExternalAngleMovement, double InternalAngleMovement, double DistanceFromRespectiveSun, double InitialOrbitalAngle, double Tilt, bool DisplayAxixPole, bool DisplayTrail, PlanetPropeties PP)
        {


            DrawingBrush ibx = new DrawingBrush();
            Type EllipisoidType = typeof(Ellipsoid);
            TubeGeometry TG;
            if (shapeType == EllipisoidType)
            {
                TG = new Ellipsoid(ModelDgroup, TranslationVector, AxisVector, TubeThickness, InternalDiameter, angleStep_horizotal, angleStep_vertical, NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, intialselfRotation, InitialgroundRotation, IsItTile, displaycommet, AxisAndCenterExternal, AxisAndCenterExternal_v2, DisplayAxixPole, DisplayTrail, ellipsisA, ellipsisB, ellipsisC, ExternalAngleMovement, InternalAngleMovement, DistanceFromRespectiveSun, InitialOrbitalAngle, Tilt, PP);
            }
            else
            {
                TG = new TubeGeometry(ModelDgroup, TranslationVector, AxisVector, TubeThickness, InternalDiameter, angleStep_horizotal, angleStep_vertical, NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, intialselfRotation, InitialgroundRotation, IsItTile, displaycommet, AxisAndCenterExternal, AxisAndCenterExternal_v2, DisplayAxixPole, DisplayTrail, ExternalAngleMovement, InternalAngleMovement, DistanceFromRespectiveSun, InitialOrbitalAngle, PP, Tilt);

            }
            ShapesList.Add(TG);
            MeshGeometry3D mesh = TG.Mesh;
            MaterialGroup MatGR = new MaterialGroup();
            MaterialGroup CircelMaterialGroup = new MaterialGroup();
            MaterialGroup AxisPoleMaterialGroup = new MaterialGroup();
            //MaterialGroup BAckMatGR = new MaterialGroup();
            if (ib is ImageBrush)
            {
                MatGR.Children.Add(new DiffuseMaterial(ib as ImageBrush));
            }
            else
            {
                MatGR.Children.Add(new DiffuseMaterial(ib as DrawingBrush));
            }
            mGeometry = new GeometryModel3D(mesh, MatGR);
            mGeometry.BackMaterial = MatGR;
            Transform3DGroup MyTransform3DGroup = new Transform3DGroup();
            MatrixTransform3D myMatrixTransform3D = new MatrixTransform3D();


            Transform3DGroup MyTransform3DGroup2 = new Transform3DGroup();
            MatrixTransform3D myMatrixTransform3D2 = new MatrixTransform3D();

            Transform3DGroup MyTransform3DGroup3 = new Transform3DGroup();
            MatrixTransform3D myMatrixTransform3D3 = new MatrixTransform3D();

            MyTransform3DGroup.Children.Add(myMatrixTransform3D);
            MyTransform3DGroup2.Children.Add(myMatrixTransform3D2);
            MyTransform3DGroup3.Children.Add(myMatrixTransform3D3);
            mGeometry.Transform = MyTransform3DGroup;
            AxisPoleMaterialGroup.Children.Add(new DiffuseMaterial(System.Windows.Media.Brushes.Green));
            CircelMaterialGroup.Children.Add(new DiffuseMaterial(System.Windows.Media.Brushes.Green));

            if (displaycommet == true)
            {
                mGeometry2 = new GeometryModel3D(TG.CircleMesh, CircelMaterialGroup);
                mGeometry2.BackMaterial = CircelMaterialGroup;
                mGeometry2.Transform = MyTransform3DGroup2;
                TG.PopulateGeometries(mGeometry, mGeometry2);
                group2.Children.Add(TG.GeometriesAndMtrices.mGeometry2);
            }
            else
            {
                mGeometry2 = new GeometryModel3D();
                TG.PopulateGeometries(mGeometry, mGeometry2);
            }

            if (DisplayAxixPole == true)
            {
                mGeometry3 = new GeometryModel3D(TG.AxisPoleMesh, AxisPoleMaterialGroup);
                mGeometry3.BackMaterial = AxisPoleMaterialGroup;
                mGeometry3.Transform = MyTransform3DGroup3;
                TG.GeometriesAndMtrices.mGeometry3 = mGeometry3;
                group3.Children.Add(TG.GeometriesAndMtrices.mGeometry3);
            }







            group.Children.Add(TG.GeometriesAndMtrices.mGeometry);


            object ffff = model;
            object zzzz = ffff;
        }
        private Model3DGroup ReturnModel3DGroup()
        {
            ModelVisual3D extraCube = new ModelVisual3D();
            Model3DGroup modelGroup = new Model3DGroup();
            GeometryModel3D model3d = new GeometryModel3D();
            extraCube.Content = modelGroup;

            viewport.Children.Add(extraCube);

            return modelGroup;
        }
        private void RecalculateTheAxis(double groundrotatingangle)
        { 
            /// the tube rotated around GroundAxis by  degrees
            /// therefore it is necessary to rotate OriginalAxis with respect to GroundAxis By groundrotatingangle
            ///  that is done by setting TheAxis to the OriginalAxis rotated by groundrotatingangle
            
        }
        private static Matrix3D returnMatrixPertransformation(double angle, Vector3D TheAxis, Vector3D TranslationVector)
        {
            TranslateTransform3D TransTran = new TranslateTransform3D();
            Transform3DGroup Tr3D = new Transform3DGroup();
            RotateTransform3D myRotateTransform3D = new RotateTransform3D();
            AxisAngleRotation3D myAxisAngleRotation3d = new AxisAngleRotation3D();
            myAxisAngleRotation3d.Axis = TheAxis;
            myAxisAngleRotation3d.Angle = angle;
            myRotateTransform3D.Rotation = myAxisAngleRotation3d;
            TransTran.OffsetX=TranslationVector.X;
            TransTran.OffsetY=TranslationVector.Y;
            TransTran.OffsetZ = TranslationVector.Z;

            Tr3D.Children.Add(TransTran);

            Tr3D.Children.Add(myRotateTransform3D);

            //mGeometry.Transform = Tr3D;
            return Tr3D.Value;
        }

        
    }
}
