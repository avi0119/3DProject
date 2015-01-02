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
    public class Point3DOriginalAndNow
    {
        protected Point3D prOriginal;
        protected Point3D prNow;
        protected Point3D prOffset;
        public Point3D Offset
        {
            get
            {
                return prOffset;
            }
            set
            {
                prOffset = value;
            }
        }
        public Point3D Now
        {
            get
            {
                return prNow;
            }
            set
            {
                
                prNow = value;
                Offset = ((Vector3D)prNow) - Original;
            }
        }
        public Point3D Original
        {
            get
            {
                return prOriginal;
            }
            set
            {
                prOriginal = value;
            }
        }
    }
    public class AxisAndCenter
    {
        protected double prCurrentAngle_self = 0;
        protected double prCurrentAngle=0;
        protected double prAngleMovemente;
        protected double prAngleMovemente_self;
        protected RotationMatrix prRotationAroundArbAxix;
        protected Vector3D prAxis;
        protected Point3D prCenter;
        protected double prTilt;
        protected double prTiltRelatedRadius;
        public double TiltRelatedRadius
        {
            get
            {
                return prTiltRelatedRadius;
            }
            set
            {
                prTiltRelatedRadius = value;
            }
        }
        public double Tilt
        {
            get
            {
                return prTilt;
            }
            set
            {
                prTilt = value;
            }
        }

        public RotationMatrix RotationAroundArbAxix
        {
            get
            {
                return prRotationAroundArbAxix;
            }
            set
            {
                prRotationAroundArbAxix = value;
            }
        }
        public Vector3D Axis
        {
            get
            {
                return prAxis;
            }
            set
            {
                prAxis = value;
            }
        }
        public Point3D Center
        {
            get
            {
                return prCenter;
            }
            set
            {
                prCenter = value;
            }
        }
        public AxisAndCenter()
        {

        }
        public double NextAngle
        {
            get
            {
                return prCurrentAngle;
            }
            set
            {
                prCurrentAngle = value;
            }
        }
        public double NextAngle_self
        {
            get
            {
                return prCurrentAngle_self;
            }
            set
            {
                prCurrentAngle_self = value;
            }
        }
        public double AngleMovemente
        {
            get
            {
                return prAngleMovemente;
            }
            set
            {
                prAngleMovemente = value;
            }
        }
        public double AngleMovemente_self
        {
            get
            {
                return prAngleMovemente_self;
            }
            set
            {
                prAngleMovemente_self = value;
            }
        }
        public AxisAndCenter(double tilt, double angleOffset, Vector3D axis, double angleOffset_self, Point3D CenterPoint = new Point3D(), double SinRotorOverCosRotor = 1,double radius=0)
        {
            RotationMatrix A = new RotationMatrix(axis, CenterPoint, SinRotorOverCosRotor);
            RotationAroundArbAxix = A;
            AngleMovemente = angleOffset;
            AngleMovemente_self = angleOffset_self;
            Center = CenterPoint;
            Tilt = tilt;
            TiltRelatedRadius = radius;
            //CalculateMatrix(AngleOfMovement);
        }
        public void CalculateMatrix (double AngleOfMovement)
        {
            RotationAroundArbAxix.Calculate3DMatrix(AngleOfMovement,0,Tilt,TiltRelatedRadius);
        }
        virtual public void CalculateMatrix()
        {

            RotationAroundArbAxix.Calculate3DMatrix(NextAngle, 0, Tilt, TiltRelatedRadius);
            NextAngle = NextAngle + AngleMovemente;
        }
        public void ApplyMAtrixToDixtionaryOfPoints(Dictionary<string, Point3DOriginalAndNow> pointsDixtionary)
        {
            foreach (string PointName  in pointsDixtionary.Keys)
            {
                Point3DOriginalAndNow CurrentOriginalAndNow = pointsDixtionary[PointName];
                CurrentOriginalAndNow.Now = ApplyToPoint3D(CurrentOriginalAndNow.Now);

            }
        }
        public Point3D ApplyToPoint3D (Point3D point)
        {


            return RotationAroundArbAxix.RotationMatrix3D.Transform(point);
          
        }


    }
    public class AxisAndCenter_v2 : AxisAndCenter
    {
        public AxisAndCenter_v2(double tilt, double angleOffset, Vector3D axis, double angleOffset_self, Point3D CenterPoint = new Point3D(), double SinRotorOverCosRotor = 1, double radius=0)//:base( angleOffset,  axis,  CenterPoint,  SinRotorOverCosRotor)
        {
            RotationMatrix A = new RotationMatrix_v2(axis, CenterPoint, SinRotorOverCosRotor);
            RotationAroundArbAxix = A;
            AngleMovemente = angleOffset;
            AngleMovemente_self = angleOffset_self;
            Center = CenterPoint;
            Axis = axis;
            Tilt = tilt;
            TiltRelatedRadius = radius;
        }
        public AxisAndCenter_v2():base()
        {

        }
        override public void CalculateMatrix()
        {
            //double tilt=0;
            RotationAroundArbAxix.Calculate3DMatrix(NextAngle, NextAngle_self, Tilt, TiltRelatedRadius);
            NextAngle = NextAngle + AngleMovemente;
            NextAngle_self=NextAngle_self+AngleMovemente_self;
        }
    }
    public class AngleVectorArray
    {
        public double Angle { get; set; }
        public Vector3D Vector { get; set; }
        public AngleVector[] array { get; set; }
    }
    public class GeometryModel3DAndMatrices
    {

        public GeometryModel3D mGeometry { get; set; }
        public GeometryModel3D mGeometry2 { get; set; }
        public GeometryModel3D mGeometry3 { get; set; }
    }
    public class AngleVector
    {
        public double Angle { get; set; }
        public Vector3D Vector { get; set; }

    }
     class TubeGeometry
     {
         protected bool prDisplayTrail;
         protected Model3DGroup prMyModel3DGroup;
         protected Vector3D prAxisToRotateAround;
         public Dictionary<string, Point3DOriginalAndNow> Important3DPoints = new Dictionary<string, Point3DOriginalAndNow>();
         // THESE POINTS ARE ROTATED WITH VERSION 2 
         public Dictionary<string, Point3DOriginalAndNow> Important3DPoints_V2 = new Dictionary<string, Point3DOriginalAndNow>();
         protected AxisAndCenter prExternalAxis_v2;
         protected AxisAndCenter prExternalAxis;
         protected  bool prDIsplayCommet;
         protected bool prDisplayPole;
         protected bool prIsItTile;
         protected MeshGeometry3D prCircleMeshMesh;
         protected MeshGeometry3D prAxisPoleMesh;
         protected MeshGeometry3D prMesh;
         protected Vector3D prAxis;
         protected Vector3D prOriginalAxis;
         protected Vector3D prGroundAxis;
         protected Vector3D prTranslationVector;
         protected Vector3D prAxisProjection;
         protected Trail3DPoints prTrail;
         protected PlanetPropeties prMyPlanetProperties;

         public PlanetPropeties MyPlanetProperties
         {
             get
             {
                 return prMyPlanetProperties;
             }
             set
             {
                 prMyPlanetProperties = value;
             }
         }

         public bool DisplayTrailFlag
         {
             get
             {
                 return prDisplayTrail;
             }
             set
             {
                 prDisplayTrail = value;
             }
         }
         public Trail3DPoints Trail
         {
             get
             {
                 return prTrail;
             }
             set
             {
                 prTrail = value;
             }
         }

         public Vector3D AxisProjection
         {
             get
             {
                 return prAxisProjection;
             }
             set
             {
                 prAxisProjection = value;
             }
         }

         //protected AxisAngleRotation3D prrotationY = new AxisAngleRotation3D();
         //protected MatrixTransform3D prmyMatrixTransform3D;
          protected double prselfrotatingangle = 0;
          protected double prgroundrotatingangle = 0;
          protected int prangleStep;
          protected double prRadius;
          protected double prTubeThicknes;
          protected GeometryModel3DAndMatrices prGeometriesAndMtrices;

          virtual public Vector3D calculateAxixProjection(Vector3D centerPointOfSocalledSun, Vector3D centerPointOfObjectInQuestion,Vector3D anAxis)
          {
              Vector3D a = anAxis / anAxis.Length;
              Vector3D v = centerPointOfObjectInQuestion - centerPointOfSocalledSun;
              v = v / v.Length;
              Vector3D N;
              //double ets2;
              if ((a - v).Length < .01)
              {
                  // the axis of the planet and the axix of the orbit are the same -- there is ni way ti tell what is the orbit
                  N = new Vector3D(a.Y, -a.X, 0);
                  //ets2 = -1;
                  
              }
              else
              {
                  double vLengthSuared = v.Length * v.Length;
                  double AdotV = Vector3D.DotProduct(a, v);
                   N = new Vector3D(v.X * AdotV - a.X * vLengthSuared, v.Y * AdotV - a.Y * vLengthSuared, v.Z * AdotV - a.Z * vLengthSuared);
                   //double test = Vector3D.DotProduct(N, a) - Vector3D.DotProduct(N, v);
                   //    ets2=test;
              }
              //double test3 = ets2;
              return N / N.Length;
          }
          //
          public Model3DGroup MyModel3DGroup
          {
              get
              {
                  return prMyModel3DGroup;
              }
              set
              {
                  prMyModel3DGroup = value;
              }
          }
          public AxisAndCenter  ExternalAxis
          {
              get
              {
                  return prExternalAxis;
              }
              set
              {
                  prExternalAxis = value;
              }
          }
          public AxisAndCenter ExternalAxis_v2
          {
              get
              {
                  return prExternalAxis_v2;
              }
              set
              {
                  prExternalAxis_v2 = value;
              }
          }

         public double RadiusOfShape
          {
              get
              {
                  return prRadius;
              }
              set
              {
                  prRadius = value;
              }
          }
         public double TUBETHICKNESS
         {
             get
             {
                 return prTubeThicknes;
             }
             set
             {
                 prTubeThicknes = value;
             }
         }
         
         public bool IsItTile
         {
             get
             {
                 return prIsItTile;
             }
             set
             {
                 prIsItTile = value;
             }
         }

         public bool DisplayCommet
         {
             get
             {
                 return prDIsplayCommet;
             }
             set
             {
                 prDIsplayCommet = value;
             }
         }
         public bool DisplayPole
         {
             get
             {
                 return prDisplayPole;
             }
             set
             {
                 prDisplayPole = value;
             }
         }

         public GeometryModel3DAndMatrices GeometriesAndMtrices
         {
             get
             {
                 return prGeometriesAndMtrices;
             }
             set
             {
                 prGeometriesAndMtrices = value;
             }
         }

         public int AngleStepPerShape
          {
              get
              {
                  return prangleStep;
              }
              set
              {
                  prangleStep = value;
              }
          }

          public MeshGeometry3D Mesh
          {
              get
              {
                  return prMesh;
              }
              set
              {
                  prMesh = value;
              }
          }
          public MeshGeometry3D CircleMesh
          {
              get
              {
                  return prCircleMeshMesh;
              }
              set
              {
                  prCircleMeshMesh = value;
              }
          }
          public MeshGeometry3D AxisPoleMesh
          {
              get
              {
                  return prAxisPoleMesh;
              }
              set
              {
                  prAxisPoleMesh = value;
              }
          }

          public double groundrotatingangle
          {
              get
              {
                  return prgroundrotatingangle;
              }
              set
              {
                  prgroundrotatingangle = value;
              }
          }
          public double selfrotatingangle
          {
              get
              {
                  return prselfrotatingangle;
              }
              set
              {
                  prselfrotatingangle = value;
              }
          }
          public Vector3D AxisToRotateAround
          {
              get
              {
                  return prAxisToRotateAround;
              }
              set
              {
                  prAxisToRotateAround = value;
              }
          }

          public Vector3D TranslationVector
          {
              get
              {
                  return prTranslationVector;
              }
              set
              {
                  prTranslationVector = value;
              }
          }
        public  Vector3D TheAxis
        {
            get
            {
                return prAxis;
            }
            set
            {
                prAxis = value;
            }
        }
        public  Vector3D OriginalAxis
        {
            get
            {
                return prOriginalAxis;
            }
            set
            {
                prOriginalAxis = value;
            }
        }
        public  Vector3D GroundAxis
        {
            get
            {
                return prGroundAxis;
            }
            set
            {
                prGroundAxis = value;
                ReturnMeshPerComet(prGroundAxis, AngleStepPerShape, RadiusOfShape);
            }
        }
        public TubeGeometry()
        {

        }
        public void PopulateGeometries(GeometryModel3D geometry1, GeometryModel3D geometry2)
        {
            this.GeometriesAndMtrices.mGeometry = geometry1;
            this.GeometriesAndMtrices.mGeometry2 = geometry2;
         }
        public TubeGeometry(Model3DGroup ModelDgroup, Vector3D translationVector_nouse, Vector3D AxisVector, double TubeThickness, double InternalDiameter, int angleStep_horizotal, int angleStep_vertical, int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, double initialSelfRotation, double intiallGroundAngle, bool isittile, bool ShoulDisplayCommet,AxisAndCenter aExternalAxis,AxisAndCenter aExternalAxis_v2,bool DisplayAxixPole,bool DisplayTrail, double ExternalAngleMovement, double InternalAngleMovemen,double DistanceFromRespectiveSun,double  InitialOrbitalAngle, PlanetPropeties PP,  double Tilt, double elilipssiA = 0, double elilipssiB = 0, double elilipssiC = 0, bool isItcomet=false)
        {
            if (GlobalVar.USEVERSION2 == true)
            {
                TubeConstructorFunction_v2(ModelDgroup, translationVector_nouse, AxisVector, TubeThickness, InternalDiameter, angleStep_horizotal, angleStep_vertical, NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, initialSelfRotation, intiallGroundAngle, isittile, ShoulDisplayCommet, aExternalAxis, aExternalAxis_v2, DisplayAxixPole, DisplayTrail, ExternalAngleMovement, InternalAngleMovemen, DistanceFromRespectiveSun, InitialOrbitalAngle, PP, Tilt, elilipssiA, elilipssiB, elilipssiC, isItcomet);
            }
            else
            {
                TubeConstructorFunction(ModelDgroup, translationVector_nouse, AxisVector, TubeThickness, InternalDiameter, angleStep_horizotal, angleStep_vertical, NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, initialSelfRotation, intiallGroundAngle, isittile, ShoulDisplayCommet, aExternalAxis, aExternalAxis_v2, DisplayAxixPole, DisplayTrail, ExternalAngleMovement, InternalAngleMovemen, DistanceFromRespectiveSun, InitialOrbitalAngle, PP, Tilt, elilipssiA, elilipssiB, elilipssiC, isItcomet);
            }
        }

        protected void TubeConstructorFunction_v2(Model3DGroup ModelDgroup, Vector3D translationVector_nouse, Vector3D AxisVector, double TubeThickness, double InternalDiameter, int angleStep_horizotal, int angleStep_vertical, int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, double initialSelfRotation, double intiallGroundAngle, bool isittile, bool ShoulDisplayCommet, AxisAndCenter aExternalAxis, AxisAndCenter aExternalAxis_v2, bool DisplayAxixPole, bool DisplayTrail, double ExternalAngleMovement, double InternalAngleMovemen, double DistanceFromRespectiveSun, double InitialOrbitalAngle, PlanetPropeties PP, double Tilt, double elilipssiA = 0, double elilipssiB = 0, double elilipssiC = 0, bool isItcomet = false)
        {
            MyPlanetProperties = PP;
            DisplayTrailFlag = DisplayTrail;
            MyModel3DGroup = ModelDgroup;
            DisplayPole = DisplayAxixPole;
            Trail = new Trail3DPoints(20, 0.0001, ModelDgroup);
            //Vector3D translationVector = CalculateTranslationVectorPerGivenCenter(AxisVector, DistanceFromRespectiveSun);
            Vector3D translationVector = (Vector3D)PP.initilaLocation / MainWindow.AstronomicalUnint;
             //translationVector = new Vector3D(); // this is for debugging the odd problem
            //CreateCubeWithinGivenModel3DGroup(MyModel3DGroup, translationVector.X, translationVector.Y, translationVector.Z,0.25);
            DisplayCommet = ShoulDisplayCommet;
            this.IsItTile = isittile;
            selfrotatingangle = initialSelfRotation;
            groundrotatingangle = intiallGroundAngle;
            GeometriesAndMtrices = new GeometryModel3DAndMatrices();
            RadiusOfShape = TubeThickness + InternalDiameter / 2.0;
            AngleStepPerShape = angleStep_horizotal;
            TranslationVector = translationVector;
            Vector3DCollection normals = new Vector3DCollection();
            Point3DCollection P3DcollectionPositions = new Point3DCollection();
            Int32Collection Triangles_int = new Int32Collection();
            PointCollection textCoord = new PointCollection();
            PointCollection textCoord_alternative = new PointCollection();

            //TG = new TubeGeometry();
            CalculateGroundAxixGivenTubePerpendicularAxis(AxisVector);

            ReturnPositionsTriangleIndicesNormals(NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, angleStep_horizotal, angleStep_vertical, TubeThickness, InternalDiameter, AxisVector, ref  P3DcollectionPositions, ref  Triangles_int, ref  normals, ref textCoord, ref  textCoord_alternative, elilipssiA, elilipssiB, elilipssiC, isItcomet);
            MeshGeometry3D mesh = new MeshGeometry3D();

            mesh.Positions = P3DcollectionPositions;
            mesh.Normals = normals;
            mesh.TriangleIndices = Triangles_int;

            mesh.TextureCoordinates = textCoord_alternative;
            Mesh = mesh;
            Important3DPoints.Add("CENTER POINT", new Point3DOriginalAndNow { Original = (Point3D)translationVector, Now = (Point3D)translationVector });
            Important3DPoints.Add("BASE AXIS", new Point3DOriginalAndNow { Original = (Point3D)translationVector, Now = (Point3D)translationVector });
            Important3DPoints.Add("END AXIS", new Point3DOriginalAndNow { Original = (Point3D)(translationVector + AxisVector), Now = (Point3D)(translationVector + AxisVector) });


            Important3DPoints_V2.Add("CENTER POINT", new Point3DOriginalAndNow { Original = (Point3D)translationVector, Now = (Point3D)translationVector });
            Important3DPoints_V2.Add("BASE AXIS", new Point3DOriginalAndNow { Original = (Point3D)translationVector, Now = (Point3D)translationVector });
            Important3DPoints_V2.Add("END AXIS", new Point3DOriginalAndNow { Original = (Point3D)(translationVector + AxisVector), Now = (Point3D)(translationVector + AxisVector) });
            TUBETHICKNESS = TubeThickness;

            //double TiltRelevantRadius = TubeThickness / 2;




            //aExternalAxis.CalculateMatrix(ExternalAngleMovement);
            ExternalAxis = aExternalAxis;
            ExternalAxis_v2 = aExternalAxis_v2;
            ExternalAxis_v2.NextAngle = InitialOrbitalAngle;
            ExternalAxis_v2.RotationAroundArbAxix.Axis = AxisToRotateAround;
            RotationMatrix_v2 aRotatiomMatrix_V2;
            aRotatiomMatrix_V2 = (RotationMatrix_v2)(ExternalAxis_v2.RotationAroundArbAxix);
            aRotatiomMatrix_V2.Important3DPoints = Important3DPoints_V2;
            AxisProjection = calculateAxixProjection((Vector3D)ExternalAxis.Center, TranslationVector, AxisVector);
            CreateAxisPoleMesh(AxisToRotateAround, TUBETHICKNESS, elilipssiA);
            //AddTranslationVecotrToAllVeritces();
        }
         private void AddTranslationVecotrToAllVeritces() 
         {
             MainWindow.globalviewport = GlobalVar.mainwindow.viewport;
             MainWindow.createGlobalCube();
             Mesh = MainWindow.globalMesh;
             var colecction=Mesh.Positions;
             //TranslationVector = new Vector3D();
             var t=colecction.Select(i=> i=i+TranslationVector);
             colecction = Mesh.Positions = new Point3DCollection( t.ToArray());

 
         }
        protected void TubeConstructorFunction(Model3DGroup ModelDgroup, Vector3D translationVector_nouse, Vector3D AxisVector, double TubeThickness, double InternalDiameter, int angleStep_horizotal, int angleStep_vertical, int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, double initialSelfRotation, double intiallGroundAngle, bool isittile, bool ShoulDisplayCommet, AxisAndCenter aExternalAxis, AxisAndCenter aExternalAxis_v2, bool DisplayAxixPole, bool DisplayTrail, double ExternalAngleMovement, double InternalAngleMovemen, double DistanceFromRespectiveSun, double InitialOrbitalAngle, PlanetPropeties PP, double Tilt, double elilipssiA = 0, double elilipssiB = 0, double elilipssiC = 0, bool isItcomet = false)
        {
            MyPlanetProperties = PP;
            DisplayTrailFlag = DisplayTrail;
            MyModel3DGroup = ModelDgroup;
            DisplayPole = DisplayAxixPole;
            Trail = new Trail3DPoints(2, 0.01, ModelDgroup);
            Vector3D translationVector = CalculateTranslationVectorPerGivenCenter(AxisVector, DistanceFromRespectiveSun);
            //CreateCubeWithinGivenModel3DGroup(MyModel3DGroup, translationVector.X, translationVector.Y, translationVector.Z,0.25);
            DisplayCommet = ShoulDisplayCommet;
            this.IsItTile = isittile;
            selfrotatingangle = initialSelfRotation;
            groundrotatingangle = intiallGroundAngle;
            GeometriesAndMtrices = new GeometryModel3DAndMatrices();
            RadiusOfShape = TubeThickness + InternalDiameter / 2.0;
            AngleStepPerShape = angleStep_horizotal;
            TranslationVector = translationVector;
            Vector3DCollection normals = new Vector3DCollection();
            Point3DCollection P3DcollectionPositions = new Point3DCollection();
            Int32Collection Triangles_int = new Int32Collection();
            PointCollection textCoord = new PointCollection();
            PointCollection textCoord_alternative = new PointCollection();

            //TG = new TubeGeometry();
            CalculateGroundAxixGivenTubePerpendicularAxis(AxisVector);

            ReturnPositionsTriangleIndicesNormals(NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, angleStep_horizotal, angleStep_vertical, TubeThickness, InternalDiameter, AxisVector, ref  P3DcollectionPositions, ref  Triangles_int, ref  normals, ref textCoord, ref  textCoord_alternative, elilipssiA, elilipssiB, elilipssiC, isItcomet);
            MeshGeometry3D mesh = new MeshGeometry3D();

            mesh.Positions = P3DcollectionPositions;
            mesh.Normals = normals;
            mesh.TriangleIndices = Triangles_int;

            mesh.TextureCoordinates = textCoord_alternative;
            Mesh = mesh;
            Important3DPoints.Add("CENTER POINT", new Point3DOriginalAndNow { Original = (Point3D)translationVector, Now = (Point3D)translationVector });
            Important3DPoints.Add("BASE AXIS", new Point3DOriginalAndNow { Original = (Point3D)translationVector, Now = (Point3D)translationVector });
            Important3DPoints.Add("END AXIS", new Point3DOriginalAndNow { Original = (Point3D)(translationVector + AxisVector), Now = (Point3D)(translationVector + AxisVector) });


            Important3DPoints_V2.Add("CENTER POINT", new Point3DOriginalAndNow { Original = (Point3D)translationVector, Now = (Point3D)translationVector });
            Important3DPoints_V2.Add("BASE AXIS", new Point3DOriginalAndNow { Original = (Point3D)translationVector, Now = (Point3D)translationVector });
            Important3DPoints_V2.Add("END AXIS", new Point3DOriginalAndNow { Original = (Point3D)(translationVector + AxisVector), Now = (Point3D)(translationVector + AxisVector) });
            TUBETHICKNESS = TubeThickness;

            //double TiltRelevantRadius = TubeThickness / 2;




            //aExternalAxis.CalculateMatrix(ExternalAngleMovement);
            ExternalAxis = aExternalAxis;
            ExternalAxis_v2 = aExternalAxis_v2;
            ExternalAxis_v2.NextAngle = InitialOrbitalAngle;
            ExternalAxis_v2.RotationAroundArbAxix.Axis = AxisToRotateAround;
            RotationMatrix_v2 aRotatiomMatrix_V2;
            aRotatiomMatrix_V2 = (RotationMatrix_v2)(ExternalAxis_v2.RotationAroundArbAxix);
            aRotatiomMatrix_V2.Important3DPoints = Important3DPoints_V2;
            AxisProjection = calculateAxixProjection((Vector3D)ExternalAxis.Center, TranslationVector, AxisVector);
            CreateAxisPoleMesh(AxisToRotateAround, TUBETHICKNESS, elilipssiA);
        }
         protected Vector3D CalculateTranslationVectorPerGivenCenter_v2(Vector3D anAxis, double LenghthOfRadius,Point3D aCenter,Plane IntersectingPlane)
        {


            return CalculateTranslationVectorPerGivenCenter_inaction_v2(anAxis, LenghthOfRadius, aCenter, IntersectingPlane);

        }
          protected Vector3D CalculateTranslationVectorPerGivenCenter_inaction_v2(Vector3D anAxis,double LenghthOfRadius,Point3D aCenter,Plane IntersectingPlane)
        {
            Plane aPlaneOfOrbit=new Plane(anAxis,aCenter);
            LineVector LV = Plane.ReturnIntersectionLine(IntersectingPlane, aPlaneOfOrbit);
            return LV.VectorDirection/LV.VectorDirection.Length*LenghthOfRadius;
        }
         //static public DrawingBrush CreateTransparemtFramedBrush()
        protected Vector3D CalculateTranslationVectorPerGivenCenter(Vector3D anAxis, double LenghthOfRadius)
        {


            return CalculateTranslationVectorPerGivenCenter_inaction( anAxis, LenghthOfRadius);

        }
        protected Vector3D CalculateTranslationVectorPerGivenCenter_inaction(Vector3D anAxis,double LenghthOfRadius)
        {
            Vector3D r1;
            Vector3D r2;
            double x; double y; double z = 0;
            x = 1; y = 1;
            Vector3D resultantInitialVector = new Vector3D(); ;

            bool resper = TryPerpendicularGivenVectorAndTwoPoints(anAxis / anAxis.Length, x, y, z, ref resultantInitialVector);
            r1 = resultantInitialVector / resultantInitialVector.Length * LenghthOfRadius;
            Vector3D b = resultantInitialVector / resultantInitialVector.Length;// new Vector3D(x / b_temp.Length, y / b_temp.Length, z / b_temp.Length);
            return r1;

            //perpendicular to axis and r1
             //r2 = Vector3D.CrossProduct(anAxis, r1);

             //r2 = r2 / r2.Length * LenghthOfRadius;

        }
         static public DrawingBrush CreateTransparemtFramedBrush()
        {
            // 
            // Create the Geometry to draw. 
            //
            GeometryGroup aRectangle = new GeometryGroup();
            aRectangle.Children.Add(new RectangleGeometry(new Rect(0, 0, 200, 200)));
            //aRectangle.Children.Add(
            //    new EllipseGeometry(new Point(50, 50), 20, 45)
            //    );

            // 
            // Create a GeometryDrawing. 
            //
            GeometryDrawing aGeometryDrawing = new GeometryDrawing();
            aGeometryDrawing.Geometry = aRectangle;

            // Paint the drawing with a gradient.
            //aGeometryDrawing.Brush =
            //    new LinearGradientBrush(
            //        Colors.Blue,
            //        Color.FromRgb(204, 204, 255),
            //        new Point(0, 0),
            //        new Point(1, 1));
            aGeometryDrawing.Brush = Brushes.Transparent;

            // Outline the drawing with a solid color.
            aGeometryDrawing.Pen = new Pen(Brushes.Blue, 10);


            DrawingGroup checkersDrawingGroup = new DrawingGroup();
            checkersDrawingGroup.Children.Add(aGeometryDrawing);
            // 
            // Use a DrawingImage and an Image control 
            // to display the drawing. 
            //
            DrawingImage geometryImage = new DrawingImage(aGeometryDrawing);

            // Freeze the DrawingImage for performance benefits.
            geometryImage.Freeze();

            DrawingBrush anImage = new DrawingBrush();
            anImage.Drawing = checkersDrawingGroup;

            return anImage;
        }
         static public DrawingBrush CreateCheckersDrawingBrush()
        {
            // Create a background recntangle
            //Rectangle chessBoard = new Rectangle();
            //chessBoard.Width = 300;
            //chessBoard.Height = 300;

            // Create a DrawingBrush
            DrawingBrush blackBrush = new DrawingBrush();
            // Create a Geometry with white background
            GeometryDrawing backgroundSquare =
                new GeometryDrawing(
                    Brushes.White,
                    null,
                    new RectangleGeometry(new Rect(0, 0, 400, 400)));
            // Create a GeometryGroup that will be added to Geometry
            GeometryGroup gGroup = new GeometryGroup();
            gGroup.Children.Add(new RectangleGeometry(new Rect(0, 0, 200, 200)));
            gGroup.Children.Add(new RectangleGeometry(new Rect(200, 200, 200, 200)));
            // Create a GeomertyDrawing
            GeometryDrawing checkers = new GeometryDrawing(new SolidColorBrush(Colors.Black), null, gGroup);

            DrawingGroup checkersDrawingGroup = new DrawingGroup();
            checkersDrawingGroup.Children.Add(backgroundSquare);
            checkersDrawingGroup.Children.Add(checkers);

            blackBrush.Drawing = checkersDrawingGroup;

            // Set Viewport and TimeMode
            blackBrush.Viewport = new Rect(0, 0, 1, 1);
            //blackBrush.TileMode = TileMode.Tile;
            return blackBrush;

            // Fill rectangle with a DrawingBrush
            //chessBoard.Fill = blackBrush;

            //LayoutRoot.Children.Add(chessBoard);
        }
         public LinearGradientBrush ReturnLinearBrush()
        {
            // Create a Rectangle
            //Rectangle blueRectangle = new Rectangle();
            //blueRectangle.Height = 100;
            //blueRectangle.Width = 200;

            // Create a linear gradient brush with five stops 
            LinearGradientBrush fiveColorLGB = new LinearGradientBrush();
            fiveColorLGB.StartPoint = new Point(0, 0);
            fiveColorLGB.EndPoint = new Point(1, 1);

            // Create and add Gradient stops
            GradientStop blueGS = new GradientStop();
            blueGS.Color = Colors.Blue;
            blueGS.Offset = 0.0;
            fiveColorLGB.GradientStops.Add(blueGS);

            GradientStop orangeGS = new GradientStop();
            orangeGS.Color = Colors.Orange;
            orangeGS.Offset = 0.25;
            fiveColorLGB.GradientStops.Add(orangeGS);

            GradientStop yellowGS = new GradientStop();
            yellowGS.Color = Colors.Brown;
            yellowGS.Offset = 0.50;
            fiveColorLGB.GradientStops.Add(yellowGS);

            GradientStop greenGS = new GradientStop();
            greenGS.Color = Colors.Green;
            greenGS.Offset = 0.75;
            fiveColorLGB.GradientStops.Add(greenGS);

            GradientStop redGS = new GradientStop();
            redGS.Color = Colors.Red;
            redGS.Offset = 1.0;
            fiveColorLGB.GradientStops.Add(redGS);
            return fiveColorLGB;

            // Set Fill property of rectangle
            //blueRectangle.Fill = fiveColorLGB;

            // Add Rectangle to the page
            //LayoutRoot.Children.Add(blueRectangle);
        }
         public void ReturnPositionsTriangleIndicesNormals(int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, int angleStep_horizontal, int angleStep_vertical, double TubeThickness, double InternalDiameter, Vector3D AxisVector, ref Point3DCollection P3DcollectionPositions, ref Int32Collection Triangles_int, ref Vector3DCollection normals, ref PointCollection textCoord, ref PointCollection textCoord_alternative, double elilipssiA, double elilipssiB, double elilipssiC,bool isItcomet )
        {

             

            TheAxis = AxisVector / AxisVector.Length;
            OriginalAxis = AxisVector / AxisVector.Length;
            AxisToRotateAround = TheAxis;
            AngleVectorArray[] VectorsArray = this.ReturnArrayOfAllPositionsInTube(NumberOfDefreesToComplete_horizontal, angleStep_horizontal, NumberOfDefreesToComplete_vertical, angleStep_vertical, TubeThickness, InternalDiameter, AxisVector,  elilipssiA,  elilipssiB,  elilipssiC);
            ReturnPositionsAndTriangleIndices(angleStep_horizontal, NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, angleStep_vertical, VectorsArray, ref P3DcollectionPositions, ref Triangles_int, ref normals, ref  textCoord, ref  textCoord_alternative, isItcomet);
        }

         virtual public void ReturnPositionsAndTriangleIndices(int angleStep_horizontal, int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, int angleStep_vertical, AngleVectorArray[] VectorsArray, ref Point3DCollection P3DcollectionPositions, ref Int32Collection Triangles_int, ref Vector3DCollection normals, ref PointCollection textCoord, ref PointCollection textCoord_alternative, bool isItcomet)
        {

            //Vector3DCollection normals = new Vector3DCollection();
            //Point3DCollection P3DcollectionPositions = new Point3DCollection(); 
            int LastPositionCounter = -1;
            List<Vector3D> Positions = new List<Vector3D>();
            List<string> Triangles = new List<string>();
            //Int32Collection Triangles_int = new Int32Collection();
            int offset = 1;
            int numberofverticalcircles = VectorsArray.Length;
            if (!((this as Ellipsoid) == null))
            {
                //++numberofverticalcircles;
            }
            else
            {
                numberofverticalcircles = VectorsArray.Length;
            }
            for (int i = 0; i < numberofverticalcircles; i = i + 1)
            {
                ReturnSetOfPositionsAndIndices(angleStep_horizontal, NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, angleStep_vertical, VectorsArray, i, offset, ref  LastPositionCounter, Positions, Triangles, Triangles_int, ref textCoord_alternative, isItcomet);
            }
            List<Point3D> PositionsPoints = (Positions.Select(a => new Point3D(a.X, a.Y, a.Z))).ToList<Point3D>();
            foreach (Vector3D v in Positions)
            {
                normals.Add(v / v.Length);
            }
            foreach (Point3D p in PositionsPoints)
            {
                P3DcollectionPositions.Add(p);
            }
            for (int counter = 0; counter < P3DcollectionPositions.Count; counter = counter + 4)
            {
                textCoord.Add(new Point(0, 0));
                textCoord.Add(new Point(0, 1));
                textCoord.Add(new Point(1, 1));
                textCoord.Add(new Point(1, 0));
            }
            double horizontalPosition = 0;
            double x = 0; double y = 0;
            for (int counter = 0; counter < NumberOfDefreesToComplete_horizontal; counter = counter + angleStep_horizontal)
            {
                double verticalPosition = 0;
                for (int internalcounter = 0; internalcounter < 360; internalcounter = internalcounter + angleStep_horizontal)
                {
                    //textCoord_alternative.Add(new Point(horizontalPosition, verticalPosition));
                    verticalPosition = verticalPosition + angleStep_horizontal / 360.0;

                }
                horizontalPosition = horizontalPosition + angleStep_horizontal / (double)NumberOfDefreesToComplete_horizontal;
            }
            //PointCollection textCoord_alternative = new PointCollection();
            //for (int counter = 0; counter < P3DcollectionPositions.Count; counter = counter + 4)
            //{
            //    textCoord_alternative.Add(new Point(0, 0));
            //    textCoord_alternative.Add(new Point(0, 1));
            //    textCoord_alternative.Add(new Point(1, 1));
            //    textCoord_alternative.Add(new Point(1, 0));
            //}
        }
         public void CalculateGroundAxixGivenTubePerpendicularAxis(Vector3D PerpendicularAxis)
        {

            double RetVectorX;
            double RetVectorY;
            double RetVectorZ;
            double PerpenX = PerpendicularAxis.X;
            double PerpenZ = PerpendicularAxis.Z;
            double PerpenY = PerpendicularAxis.Y;
            if (!(PerpenZ == 0))
            {

                RetVectorX = 1;
                RetVectorY = 0;//  we are on the ground plane

                double dotProductMultiplication_knowPart = PerpenX * RetVectorX + RetVectorY * RetVectorY;
                RetVectorZ = -dotProductMultiplication_knowPart / PerpenZ;
            }
            else if (!(PerpenY == 0))
            {
                RetVectorX = 1;
                RetVectorZ = 0;//  we are on the ground plane

                double dotProductMultiplication_knowPart = PerpenX * RetVectorX + RetVectorZ * RetVectorZ;
                RetVectorY = -dotProductMultiplication_knowPart / PerpenY;
            }
            else
            {
                RetVectorY = 1;
                RetVectorZ = 0;//  we are on the ground plane

                double dotProductMultiplication_knowPart = PerpenY * RetVectorY + RetVectorZ * RetVectorZ;
                RetVectorX = -dotProductMultiplication_knowPart / PerpenX;
            }
            Vector3D ret = new Vector3D(RetVectorX, RetVectorY, RetVectorZ);
            GroundAxis= ret / ret.Length;
        }
         virtual public void ReturnSetOfPositionsAndIndices(int angleStep_horizontal, int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, int angleStep_vertical, AngleVectorArray[] VectorsArray, int CurrentMainOrdinalAngle, int offset, ref int LastPositionCounter, List<Vector3D> Positions, List<string> Triangles, Int32Collection Triangles_int, ref PointCollection textCoord_alternative,bool isItcomet)// offset can be 1 or minus 1
        {
            //PointCollection textCoord_alternative = new PointCollection();
            int z=1;
            if (CurrentMainOrdinalAngle == 74)
            {
                 z = 0;
            }
             int rr=z;
            //if (((this as Ellipsoid) == null))
            //{
             if (isItcomet == true)
             {
                 if ((NumberOfDefreesToComplete_horizontal < 360) && (((CurrentMainOrdinalAngle + offset) * angleStep_horizontal) == NumberOfDefreesToComplete_horizontal) && NumberOfDefreesToComplete_horizontal == 360)
                 {
                     return;
                 }
                 if ((((CurrentMainOrdinalAngle + offset) * angleStep_horizontal) > NumberOfDefreesToComplete_horizontal))
                 {
                     return;

                 }
             }
             else
             {
                 if ((NumberOfDefreesToComplete_horizontal < 360) && (((CurrentMainOrdinalAngle + offset) * angleStep_horizontal) == NumberOfDefreesToComplete_horizontal))
                 {
                     return;
                 }
                 if ((((CurrentMainOrdinalAngle + offset) * angleStep_horizontal) > NumberOfDefreesToComplete_horizontal))
                 {
                     return;

                 }
             }
            //}
            AngleVectorArray SecondCircle;
            AngleVectorArray FirstCircle = VectorsArray[CurrentMainOrdinalAngle];
            if ((DisplayCommet == true ) ||(DisplayPole == true))
            {
                if (((((CurrentMainOrdinalAngle + offset) * angleStep_horizontal == NumberOfDefreesToComplete_horizontal) && (NumberOfDefreesToComplete_horizontal == 360)) ? 0 : (CurrentMainOrdinalAngle + offset)) == VectorsArray.Length)
                {
                    return;
                }
                else
                {
                    SecondCircle = VectorsArray[(((CurrentMainOrdinalAngle + offset) * angleStep_horizontal == NumberOfDefreesToComplete_horizontal) && (NumberOfDefreesToComplete_horizontal == 360)) ? 0 : (CurrentMainOrdinalAngle + offset)];
                }
            }
            else
            {
                if (!((this as Ellipsoid) == null))
                {
                    //++numberofverticalcircles;
                    if ((VectorsArray.Length>=(CurrentMainOrdinalAngle+1) ))
                    {
                        SecondCircle = VectorsArray[CurrentMainOrdinalAngle + 1];
                    }  
                    else
                    {
                        //return;
                        SecondCircle = VectorsArray[0];
                    }
                }
                else
                {
                    SecondCircle = VectorsArray[(((CurrentMainOrdinalAngle + offset) * angleStep_horizontal == NumberOfDefreesToComplete_horizontal) && (true)) ? 0 : (CurrentMainOrdinalAngle + offset)];
                }
            }
            bool flagOfend_horizontal = (CurrentMainOrdinalAngle + offset) * angleStep_horizontal == NumberOfDefreesToComplete_horizontal;
            AngleVector Firsti = new AngleVector();
            AngleVector Firstiplus1 = new AngleVector();
            AngleVector Secondi = new AngleVector();
            AngleVector Secondiplus1 = new AngleVector();
            bool continueflag;
            double marginfraction = 0;
            //double VerticalStepPerTextCoordinates = (1.0-marginfraction) / FirstCircle.array.Length;
            double VerticalStepPerTextCoordinates = (1.0 - marginfraction) / NumberOfDefreesToComplete_vertical / angleStep_vertical;
            double HorizontalStepPerTextCoordinates = (1.0 - marginfraction) / (NumberOfDefreesToComplete_horizontal / angleStep_horizontal);
            for (int i = 0; i < FirstCircle.array.Length; i = i + 1)
            {
                /*
                 -----------------------
                 ---- Positions Array---
                 -----------------------
                  first i=0 position
                  second i=1 position
                  second i+1=2 position
                  first i+1= 3 position
                  
                 */


                if (i == FirstCircle.array.Length - 1)
                {
                    Firsti = (FirstCircle.array)[i];

                    if (Firsti.Angle == (360 - angleStep_vertical))
                    {
                        Secondi = (SecondCircle.array)[i];
                        Firstiplus1 = (FirstCircle.array)[0];
                        Secondiplus1 = (SecondCircle.array)[0];
                        continueflag = true;
                    }
                    else
                    {
                        continueflag = true;
                    }
                }
                else
                {
                    Firsti = (FirstCircle.array)[i];
                    Firstiplus1 = (FirstCircle.array)[i + 1];
                    Secondi = (SecondCircle.array)[i];
                    Secondiplus1 = (SecondCircle.array)[i + 1];
                    continueflag = true;
                }
                if (continueflag == true)
                {
                    Positions.Add(Firsti.Vector);
                    Positions.Add(Secondi.Vector);
                    Positions.Add(Secondiplus1.Vector);
                    Positions.Add(Firstiplus1.Vector);

                    // first triangle positions
                    /// first i, second i, second i+1
                    Triangles.Add((LastPositionCounter + 1).ToString() + "," + (LastPositionCounter + 2).ToString() + "," + (LastPositionCounter + 3).ToString());
                    Triangles_int.Add(LastPositionCounter + 1);
                    Triangles_int.Add(LastPositionCounter + 2);
                    Triangles_int.Add(LastPositionCounter + 3);

                    //  second triangle positions
                    /// first i,  second i+1,first  i+1,
                    Triangles.Add((LastPositionCounter + 1).ToString() + "," + (LastPositionCounter + 3).ToString() + "," + (LastPositionCounter + 4).ToString());

                    Triangles_int.Add(LastPositionCounter + 1);
                    Triangles_int.Add(LastPositionCounter + 3);
                    Triangles_int.Add(LastPositionCounter + 4);
                    double x = marginfraction / 2.0 + (double)CurrentMainOrdinalAngle * 1.0 * HorizontalStepPerTextCoordinates;
                    double y = marginfraction / 2.0 + (double)i * 1.0 * VerticalStepPerTextCoordinates;

                    //  commented out
                    bool flagOfend_vertical = (Firsti.Angle == (360 - angleStep_vertical));
                    if (flagOfend_vertical == true)
                    {
                        int yyy = 0;
                    }
                    int lastcont = textCoord_alternative.Count;
                    if (this.IsItTile == false)
                    {
                        textCoord_alternative.Add(new Point(x, y));
                        textCoord_alternative.Add(new Point(flagOfend_horizontal == true ? 1 : x + HorizontalStepPerTextCoordinates, y));
                        textCoord_alternative.Add(new Point(flagOfend_horizontal == true ? 1 : x + HorizontalStepPerTextCoordinates, flagOfend_vertical == true ? 1 : y + VerticalStepPerTextCoordinates));
                        textCoord_alternative.Add(new Point(x, flagOfend_vertical == true ? 1 : y + VerticalStepPerTextCoordinates));
                    }
                    else
                    {
                        textCoord_alternative.Add(new Point(0, 0));
                        textCoord_alternative.Add(new Point(1, 0));
                        textCoord_alternative.Add(new Point(1, 1));
                        textCoord_alternative.Add(new Point(0, 1));
                    }



                    //textCoord.Add(new Point(0, 0));
                    //textCoord.Add(new Point(0, 1));
                    //textCoord.Add(new Point(1, 1));
                    //textCoord.Add(new Point(1, 0));

                    LastPositionCounter = LastPositionCounter + 4;
                    if (i == 0)
                    {
                        textCoord_alternative[lastcont] = textCoord_alternative[lastcont + 1];
                    }
                }

            }

            //object[] Result = new object[1];
            //return Result;
            return;
        }

         protected void ReturnMeshPerComet(Vector3D AnAxis, int angleStep,  double Radius)
         {
             //double imaginaryAngle = angleStep+angleStep;
             double pi = Math.PI;
             double imaginaryAngle = angleStep ;
             int NumberOfDefreesToComplete = 360;
             Vector3D InternalTranslatorVector = AnAxis / AnAxis.Length *( 49*pi*Radius * imaginaryAngle/360);
             AngleVector[] arr1 = ReturnArrayOfCirclePerGroundAxis(AnAxis,  angleStep,  NumberOfDefreesToComplete,  Radius,  InternalTranslatorVector);
             AngleVector[] arr2 = ReturnArrayOfCirclePerGroundAxis(AnAxis,  angleStep,  NumberOfDefreesToComplete,  Radius,  -InternalTranslatorVector);
             Vector3D currFinalFinal=new Vector3D();


             AngleVectorArray temp;

             AngleVectorArray[] final = null;
             final = new AngleVectorArray[2];

             temp = new AngleVectorArray { Angle = 0, Vector = currFinalFinal, array = arr1 };
             final[0] = temp;
              temp = new AngleVectorArray { Angle = imaginaryAngle, Vector = currFinalFinal, array = arr2 };
             final[1] = temp;
             Point3DCollection P3DcollectionPositions=new Point3DCollection();
             Int32Collection Triangles_int=new Int32Collection();
             Vector3DCollection normals=new Vector3DCollection();
             PointCollection textCoord=new PointCollection();
             PointCollection textCoord_alternative=new PointCollection() ;

             ReturnPositionsAndTriangleIndices((int)angleStep, (int)imaginaryAngle, 360, (int)angleStep, final, ref P3DcollectionPositions, ref Triangles_int, ref normals, ref  textCoord, ref  textCoord_alternative, true );

             MeshGeometry3D circleMesh = new MeshGeometry3D();
             circleMesh.Positions = P3DcollectionPositions;
             circleMesh.Normals = normals;
             circleMesh.TriangleIndices = Triangles_int;

             circleMesh.TextureCoordinates = textCoord_alternative;
             CircleMesh = circleMesh;
             //Cube mycube = new Cube();
             //CircleMesh = mycube.Mesh;


         }


         protected virtual void CreateAxisPoleMesh(Vector3D aAxis, double shapesReadius, double aelilipssiA)
         {
             Vector3D AnAxis = aAxis;
             Point3D shapesCenter= new Point3D(0,0,0);
             //double shapesReadius = TUBETHICKNESS;
             MeshGeometry3D ret = ReturnMeshGeometry3DPerAxisOfSelfRotation( AnAxis,  shapesCenter,  shapesReadius);
             AxisPoleMesh = ret;
         }
         protected MeshGeometry3D ReturnMeshGeometry3DPerAxisOfSelfRotation(Vector3D AnAxis, Point3D shapesCenter, double shapesReadius)
         {
             double Radius = 1;
             int angleStep = 120;
             Vector3D centerVector =(Vector3D) shapesCenter;
             //double imaginaryAngle = angleStep+angleStep;
             double pi = Math.PI;
             double imaginaryAngle = angleStep;
             int NumberOfDefreesToComplete = 360;
             Vector3D InternalTranslatorVector_beginningOfPole = centerVector;
             Vector3D InternalTranslatorVector_endOfPole = centerVector + shapesReadius * 2 * AnAxis;
             AngleVector[] arr1 = ReturnArrayOfCirclePerGroundAxis(AnAxis, angleStep, NumberOfDefreesToComplete, Radius, InternalTranslatorVector_beginningOfPole);
             AngleVector[] arr2 = ReturnArrayOfCirclePerGroundAxis(AnAxis, angleStep, NumberOfDefreesToComplete, Radius, InternalTranslatorVector_endOfPole);
             Vector3D currFinalFinal = new Vector3D();


             AngleVectorArray temp;

             AngleVectorArray[] final = null;
             final = new AngleVectorArray[2];

             temp = new AngleVectorArray { Angle = 0, Vector = currFinalFinal, array = arr1 };
             final[0] = temp;
             temp = new AngleVectorArray { Angle = imaginaryAngle, Vector = currFinalFinal, array = arr2 };
             final[1] = temp;
             Point3DCollection P3DcollectionPositions = new Point3DCollection();
             Int32Collection Triangles_int = new Int32Collection();
             Vector3DCollection normals = new Vector3DCollection();
             PointCollection textCoord = new PointCollection();
             PointCollection textCoord_alternative = new PointCollection();

             ReturnPositionsAndTriangleIndices((int)angleStep, (int)imaginaryAngle, 360, (int)angleStep, final, ref P3DcollectionPositions, ref Triangles_int, ref normals, ref  textCoord, ref  textCoord_alternative, true);

             MeshGeometry3D axisPoleMesh = new MeshGeometry3D();
             axisPoleMesh.Positions = P3DcollectionPositions;
             axisPoleMesh.Normals = normals;
             axisPoleMesh.TriangleIndices = Triangles_int;

             axisPoleMesh.TextureCoordinates = textCoord_alternative;
             return  axisPoleMesh;
             //Cube mycube = new Cube();
             //CircleMesh = mycube.Mesh;


         }

         protected AngleVector[] ReturnArrayOfCirclePerGroundAxis( Vector3D AnAxis, int angleStep, int NumberOfDefreesToComplete_vertical, double radiusPerPerpemdicularCircle, Vector3D InternalTranslatorVector)
         {
             double radiusPerImaginaryCircle=0;
             double x; double y; double z = 1;
             x = 1; y = 1;
             Vector3D resultantInitialVector = new Vector3D();;

             bool resper = TryPerpendicularGivenVectorAndTwoPoints(AnAxis, x, y, z, ref resultantInitialVector);
             Vector3D a_normalized = resultantInitialVector / resultantInitialVector.Length;
             AngleVector[] ResultArray = null;
             Vector3D b = resultantInitialVector / resultantInitialVector.Length;// new Vector3D(x / b_temp.Length, y / b_temp.Length, z / b_temp.Length);


             //perpendicular to a and b
             Vector3D c_temp = Vector3D.CrossProduct(AnAxis, b);
             //Vector3D d = new Vector3D(c_temp.X / c_temp.Length, c_temp.Y / c_temp.Length, c_temp.Z / c_temp.Length);
             Vector3D PerpendicularToThatWhichIsPerpendicularForAxis = new Vector3D(c_temp.X, c_temp.Y, c_temp.Z) / c_temp.Length;

             Vector3D PerpendicularToThatWhichIsPerpendicularForAxis_final = PerpendicularToThatWhichIsPerpendicularForAxis * radiusPerImaginaryCircle;
             Return360pointsGivenTwoPerPendicularVectors(NumberOfDefreesToComplete_vertical, angleStep, a_normalized, PerpendicularToThatWhichIsPerpendicularForAxis, ref ResultArray, radiusPerPerpemdicularCircle, InternalTranslatorVector);
             return ResultArray;
         }
         virtual public AngleVectorArray[] ReturnArrayOfAllPositionsInTube(int NumberOfDefreesToComplete_horizontal, int angleStep_horizontal, int NumberOfDefreesToComplete_vertical, int angleStep_vertical, double TubeThickness, double InternalDiameter, Vector3D AxisVector, double elilipssiA=0, double elilipssiB=0, double elilipssiC=0)
        {


            double radiusPerImaginaryCircle = InternalDiameter / 2 + TubeThickness / 2;
            double radiusPerPerpemdicularCircle = TubeThickness / 2;
            Vector3D b_temp;
            b_temp = new Vector3D();

            Vector3D a_normalized;

            a_normalized = AxisVector / AxisVector.Length;
            double x; double y; double z = 0;
            x = 1; y = 1;


            bool resper = TryPerpendicularGivenVectorAndTwoPoints(AxisVector, x, y, z, ref b_temp);


            double anglebetAandB = FindAngleBtweemTwoVectors(AxisVector, b_temp);


            // perpendicular to a passing through x=1 y=0
            Vector3D b = b_temp / b_temp.Length;// new Vector3D(x / b_temp.Length, y / b_temp.Length, z / b_temp.Length);


            //perpendicular to a and b
            Vector3D c_temp = Vector3D.CrossProduct(AxisVector, b);
            //Vector3D d = new Vector3D(c_temp.X / c_temp.Length, c_temp.Y / c_temp.Length, c_temp.Z / c_temp.Length);
            Vector3D c = new Vector3D(c_temp.X, c_temp.Y, c_temp.Z) / c_temp.Length;
            //double dp = Vector3D.DotProduct(a,b);
            double pi = Math.PI;
            double currentRadian;
            Vector3D First;
            Vector3D Second;
            double Cosine;
            Double Sine;
            Vector3D Final;
            Double dotProductTest;


            AngleVector[] arr = null;
            AngleVectorArray[] final = null;
            final = new AngleVectorArray[NumberOfDefreesToComplete_horizontal / angleStep_horizontal];
            int j = 0;
            Vector3D currFinalFinal;
            for (double i = 0; i < NumberOfDefreesToComplete_horizontal; i = i + angleStep_horizontal)
            {
                currentRadian = i / 180 * pi;
                Cosine = Math.Cos(currentRadian);
                Sine = Math.Sin(currentRadian);
                First = b * Cosine;
                Second = c * Sine;
                Final = First + Second;

                double AngleInDegrees = FindAngleBtweemTwoVectors(b, Final);
                dotProductTest = Vector3D.DotProduct(AxisVector, Final);
                double anglebetAandFinal = FindAngleBtweemTwoVectors(a_normalized, Final);
                currFinalFinal = radiusPerImaginaryCircle * Final;
                Return360pointsGivenTwoPerPendicularVectors(NumberOfDefreesToComplete_vertical, angleStep_vertical, a_normalized, Final, ref arr, radiusPerPerpemdicularCircle, currFinalFinal);
                AngleVectorArray temp = new AngleVectorArray { Angle = i, Vector = currFinalFinal, array = arr };
                final[j] = temp;
                j++;
            }
            return final;

        }
         protected bool TryPerpendicularGivenVectorAndTwoPoints(Vector3D VectorInquestion, double x, double y, double z, ref Vector3D ResultVector)
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
         protected double FindAngleBtweemTwoVectors(Vector3D a, Vector3D b)
        {
            double DotProd;
            double CosineOfAngleBetween;
            double pi = Math.PI;
            double AngleBetween;
            double AngleInDegrees;

            DotProd = Vector3D.DotProduct(b, a);
            CosineOfAngleBetween = DotProd / (b.Length * a.Length);
            AngleBetween = Math.Acos(CosineOfAngleBetween);
            AngleInDegrees = AngleBetween / pi * 180;
            return AngleInDegrees;
        }
         virtual public void Return360pointsGivenTwoPerPendicularVectors(int NumberOfDefreesToComplete_vertical, int angleStep_vertical, Vector3D a, Vector3D ResultVector, ref AngleVector[] arr, double radiusPerPerpemdicularCircle, Vector3D CenterPointVector,bool IsItForShapeItself=true)
        {
            arr = new AngleVector[NumberOfDefreesToComplete_vertical / angleStep_vertical];
            double pi = Math.PI;
            double currentRadian;
            Vector3D First;
            Vector3D Second;
            double Cosine;
            Double Sine;
            Vector3D Final;
            //Double dotProductTest;
            //double DotBFinal;
            //double CosineOfAngleBetween;
            //double AngleBetween;
            double AngleInDegrees_test;
            int j = 0;
            Vector3D currFinalFinal;
            for (double i = 0; i < NumberOfDefreesToComplete_vertical; i = i + angleStep_vertical)
            {

                currentRadian = i / 180 * pi;
                Cosine = Math.Cos(currentRadian);
                Sine = Math.Sin(currentRadian);
                First = ResultVector * Cosine;
                Second = a * Sine;
                Final = First + Second;
                //DotBFinal = Vector3D.DotProduct(ResultVector, Final);
                //CosineOfAngleBetween = DotBFinal / (ResultVector.Length * Final.Length);
                //AngleBetween = Math.Acos(CosineOfAngleBetween);
                AngleInDegrees_test = FindAngleBtweemTwoVectors(ResultVector, Final);
                currFinalFinal = Final * radiusPerPerpemdicularCircle;
                currFinalFinal = currFinalFinal + CenterPointVector;
                AngleVector temp = new AngleVector { Angle = i, Vector = currFinalFinal + TranslationVector };
                arr[j] = temp;
                j++;
                //dotProductTest = Vector3D.DotProduct(a, Final);
            }
        }

    }

}
