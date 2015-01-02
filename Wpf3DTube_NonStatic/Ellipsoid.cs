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
    class Ellipsoid : TubeGeometry
    {
        public Ellipsoid(Model3DGroup ModelDgroup,Vector3D translationVector, Vector3D AxisVector, double TubeThickness, double InternalDiameter, int angleStep_horizotal, int angleStep_vertical, int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, double initialSelfRotation, double intiallGroundAngle, bool isittile, bool ShoulDisplayCommet, AxisAndCenter ExternalAxis, AxisAndCenter ExternalAxis_v2,bool DisplayAxixPole,bool DisplayTrail, double elilipssiA, double elilipssiB, double elilipssiC, double ExternalAngleMovement, double InternalAngleMovement, double DistanceFromRespectiveSun, double InitialOrbitalAngle, double Tilt,PlanetPropeties PP)
            : base(ModelDgroup,translationVector, AxisVector, TubeThickness, InternalDiameter, angleStep_horizotal, angleStep_vertical, NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, initialSelfRotation, intiallGroundAngle, isittile, ShoulDisplayCommet, ExternalAxis, ExternalAxis_v2,DisplayAxixPole,DisplayTrail, ExternalAngleMovement, InternalAngleMovement, DistanceFromRespectiveSun, InitialOrbitalAngle,PP, Tilt, elilipssiA, elilipssiB, elilipssiC)
        {
            int z = 0;
               int r=z;
        }
        protected override void CreateAxisPoleMesh(Vector3D aAxis, double shapesReadius, double aelilipssiA)
        {
            shapesReadius = aelilipssiA;
            Vector3D AnAxis = aAxis;
            Point3D shapesCenter = new Point3D(0,0,0);
            //double shapesReadius = TUBETHICKNESS;
            MeshGeometry3D ret = ReturnMeshGeometry3DPerAxisOfSelfRotation(AnAxis, shapesCenter, shapesReadius);
            AxisPoleMesh = ret;
        }
        override public AngleVectorArray[] ReturnArrayOfAllPositionsInTube(int NumberOfDefreesToComplete_horizontal, int angleStep_horizontal, int NumberOfDefreesToComplete_vertical, int angleStep_vertical, double TubeThickness, double InternalDiameter, Vector3D AxisVector, double elilipssiA, double elilipssiB, double elilipssiC)
        {

            //NumberOfDefreesToComplete_horizontal = NumberOfDefreesToComplete_horizontal / 2;
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
            Vector3D A = a_normalized*elilipssiA; //original axix on whcih the ellipsoid spins
            Vector3D B = b * elilipssiB; /// one axis perpendiculat to A and C 
            Vector3D C = c * elilipssiC; /// second axis perpendiculat to A and B 

            AngleVector[] arr = null;
            AngleVectorArray[] final = null;
            final = new AngleVectorArray[NumberOfDefreesToComplete_horizontal / angleStep_horizontal+1];
            int j = 0;
            Vector3D currFinalFinal;
            for (double i = 0; i <= NumberOfDefreesToComplete_horizontal; i = i + angleStep_horizontal)
            {
                currentRadian = i / 180 * pi;
                Cosine = Math.Cos(currentRadian);
                Sine = Math.Sin(currentRadian);
                Vector3D VectorInDirectionOfC = C * Sine;
                First = B * Cosine;
                Second = VectorInDirectionOfC;
                
                Final = First + Second;
                double LengthOfFirstInternalEllipssisFirstAxis = Second.Length;
                Vector3D VectorInDirectionofA = A * Sine;
                double LengthOfSecondInternalEllipssisFirstAxis = VectorInDirectionofA.Length;
                // imaginary ellipsis need to be drawn and translated by 'First' 
                // one axix as the direction of A and its ;ength is  LengthOfSecondInternalEllipssisFirstAxis
                // the other is in the direction of C and its length is LengthOfFirstInternalEllipssisFirstAxis
                //  LengthOfFirstInternalEllipssisFirstAxis and LengthOfSecondInternalEllipssisFirstAxis
                // are perpendicular and are the ellipsis so-called radiuses which need to be used withoin a
                // an angle range of 360 degrees

                //double AngleInDegrees = FindAngleBtweemTwoVectors(b, Final);
                //dotProductTest = Vector3D.DotProduct(AxisVector, Final);
                //double anglebetAandFinal = FindAngleBtweemTwoVectors(a_normalized, Final);
                currFinalFinal = First;
                Return360pointsGivenTwoPerPendicularVectors(NumberOfDefreesToComplete_vertical, angleStep_vertical, VectorInDirectionofA, VectorInDirectionOfC, ref arr, radiusPerPerpemdicularCircle, currFinalFinal,true);
                AngleVectorArray temp = new AngleVectorArray { Angle = i, Vector = currFinalFinal, array = arr };
                final[j] = temp;
                j++;
            }
            TheAxis = B / B.Length;
            OriginalAxis = B / B.Length;
            return final;

        }
        override public void Return360pointsGivenTwoPerPendicularVectors(int NumberOfDefreesToComplete_vertical, int angleStep_vertical, Vector3D A, Vector3D C, ref AngleVector[] arr, double radiusPerPerpemdicularCircle, Vector3D CenterPointVector,bool IsItForShapeItself=false)
        {
            if (IsItForShapeItself == false)
            {
                base.Return360pointsGivenTwoPerPendicularVectors(NumberOfDefreesToComplete_vertical,  angleStep_vertical,  A,  C, ref arr,  radiusPerPerpemdicularCircle,  CenterPointVector,true);
                return;
            }
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
                First = C * Cosine;
                Second = A * Sine;
                Final = First + Second;
                //DotBFinal = Vector3D.DotProduct(ResultVector, Final);
                //CosineOfAngleBetween = DotBFinal / (ResultVector.Length * Final.Length);
                //AngleBetween = Math.Acos(CosineOfAngleBetween);
                AngleInDegrees_test = FindAngleBtweemTwoVectors(C, Final);
                currFinalFinal = Final * 1;
                //currFinalFinal = Final * radiusPerPerpemdicularCircle;
                currFinalFinal = currFinalFinal + CenterPointVector;
                AngleVector temp = new AngleVector { Angle = i, Vector = currFinalFinal + TranslationVector };
                arr[j] = temp;
                j++;
                //dotProductTest = Vector3D.DotProduct(a, Final);
            }
        }

        override public void ReturnPositionsAndTriangleIndices(int angleStep_horizontal, int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, int angleStep_vertical, AngleVectorArray[] VectorsArray, ref Point3DCollection P3DcollectionPositions, ref Int32Collection Triangles_int, ref Vector3DCollection normals, ref PointCollection textCoord, ref PointCollection textCoord_alternative, bool isItcomet)
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
        override  public void ReturnSetOfPositionsAndIndices(int angleStep_horizontal, int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, int angleStep_vertical, AngleVectorArray[] VectorsArray, int CurrentMainOrdinalAngle, int offset, ref int LastPositionCounter, List<Vector3D> Positions, List<string> Triangles, Int32Collection Triangles_int, ref PointCollection textCoord_alternative,bool isItcomet)// offset can be 1 or minus 1
        {
            //PointCollection textCoord_alternative = new PointCollection();
            if (CurrentMainOrdinalAngle == 3)
            {
                int z = 0;
            }
            //if (((this as Ellipsoid) == null))
            //{
            if (isItcomet == true)
            {
                if ((NumberOfDefreesToComplete_horizontal < 180) && (((CurrentMainOrdinalAngle + offset) * angleStep_horizontal) == NumberOfDefreesToComplete_horizontal) && NumberOfDefreesToComplete_horizontal == 180)
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
            if (isItcomet == true)
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
                //if (!((this as Ellipsoid) == null))
                //{
                    //++numberofverticalcircles;
                //if ((VectorsArray.Length > (CurrentMainOrdinalAngle + 1)) || ((((CurrentMainOrdinalAngle + offset) * angleStep_horizontal) < 180)))
                    if ( (VectorsArray.Length > (CurrentMainOrdinalAngle + 1)) || false )
                    {
                        SecondCircle = VectorsArray[CurrentMainOrdinalAngle + 1];
                    }
                    else
                    {
                        //return;
                        if ((CurrentMainOrdinalAngle * angleStep_horizontal) < 180) return;
                        SecondCircle = VectorsArray[0];
                    }
                //}
                //else
                //{
                //    SecondCircle = VectorsArray[(((CurrentMainOrdinalAngle + offset) * angleStep_horizontal == NumberOfDefreesToComplete_horizontal) && (true)) ? 0 : (CurrentMainOrdinalAngle + offset)];
                //}
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

    }

}
