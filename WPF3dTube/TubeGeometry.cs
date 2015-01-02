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

namespace WPF3dTube
{
    public class AngleVectorArray
    {
        public double Angle { get; set; }
        public Vector3D Vector { get; set; }
        public AngleVector[] array { get; set; }
    }
    public class AngleVector
    {
        public double Angle { get; set; }
        public Vector3D Vector { get; set; }

    }
    static class TubeGeometry
    {
        static private Vector3D prAxis;
        static private Vector3D prOriginalAxis;
        static private Vector3D prGroundAxis;

        public static Vector3D TheAxis
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
        public static Vector3D OriginalAxis
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
        public static Vector3D GroundAxis
        {
            get
            {
                return prGroundAxis;
            }
            set
            {
                prGroundAxis = value;
            }
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
            aGeometryDrawing.Brush =Brushes.Transparent;

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
        static public LinearGradientBrush ReturnLinearBrush()
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
        static public void ReturnPositionsTriangleIndicesNormals(int NumberOfDefreesToComplete_horizontal,int  NumberOfDefreesToComplete_vertical, int angleStep_horizontal,int angleStep_vertical, double TubeThickness, double InternalDiameter, Vector3D AxisVector, ref Point3DCollection P3DcollectionPositions, ref Int32Collection Triangles_int, ref Vector3DCollection normals, ref PointCollection textCoord, ref PointCollection textCoord_alternative)
        {



            TheAxis = AxisVector / AxisVector.Length;
            OriginalAxis = AxisVector / AxisVector.Length;
            AngleVectorArray[] VectorsArray = ReturnArrayOfAllPositionsInTube(NumberOfDefreesToComplete_horizontal, angleStep_horizontal, NumberOfDefreesToComplete_vertical, angleStep_vertical,TubeThickness, InternalDiameter, AxisVector);
            ReturnPositionsAndTriangleIndices(angleStep_horizontal,NumberOfDefreesToComplete_horizontal, NumberOfDefreesToComplete_vertical, angleStep_vertical, VectorsArray, ref P3DcollectionPositions, ref Triangles_int, ref normals,ref  textCoord,ref  textCoord_alternative);
        }
        static public void ReturnPositionsAndTriangleIndices(int angleStep_horizontal, int NumberOfDefreesToComplete_horizontal,int  NumberOfDefreesToComplete_vertical,int angleStep_vertical, AngleVectorArray[] VectorsArray, ref Point3DCollection P3DcollectionPositions, ref Int32Collection Triangles_int, ref Vector3DCollection normals, ref PointCollection textCoord, ref PointCollection textCoord_alternative)
        {
            //Vector3DCollection normals = new Vector3DCollection();
            //Point3DCollection P3DcollectionPositions = new Point3DCollection(); 
            int LastPositionCounter = -1;
            List<Vector3D> Positions = new List<Vector3D>();
            List<string> Triangles = new List<string>();
            //Int32Collection Triangles_int = new Int32Collection();
            int offset = 1;
            for (int i = 0; i < VectorsArray.Length; i = i + 1)
            {
                ReturnSetOfPositionsAndIndices(angleStep_horizontal,NumberOfDefreesToComplete_horizontal,  NumberOfDefreesToComplete_vertical, angleStep_vertical, VectorsArray, i, offset, ref  LastPositionCounter, Positions, Triangles, Triangles_int,ref textCoord_alternative);
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
                textCoord.Add(new Point(0,0));
                textCoord.Add(new Point(0, 1));
                textCoord.Add(new Point(1, 1));
                textCoord.Add(new Point(1, 0));
            }
            double horizontalPosition=0;
            double x = 0; double y = 0;
            for (int counter = 0; counter < NumberOfDefreesToComplete_horizontal; counter = counter + angleStep_horizontal)
            {
                double verticalPosition = 0;
                for (int internalcounter = 0; internalcounter < 360; internalcounter = internalcounter + angleStep_horizontal)
                {
                    //textCoord_alternative.Add(new Point(horizontalPosition, verticalPosition));
                    verticalPosition = verticalPosition + angleStep_horizontal / 360.0;

                }
                horizontalPosition = horizontalPosition + angleStep_horizontal /( double)NumberOfDefreesToComplete_horizontal;
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
        static public Vector3D FindGroundAxixGivenTubePerpendicularAxis(Vector3D PerpendicularAxis)
        {

            double RetVectorX ;
            double RetVectorY ; 
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
            return  ret/ret.Length;
        }
        static public void ReturnSetOfPositionsAndIndices(int angleStep_horizontal, int NumberOfDefreesToComplete_horizontal, int NumberOfDefreesToComplete_vertical, int angleStep_vertical, AngleVectorArray[] VectorsArray, int CurrentMainOrdinalAngle, int offset, ref int LastPositionCounter, List<Vector3D> Positions, List<string> Triangles, Int32Collection Triangles_int, ref PointCollection textCoord_alternative)// offset can be 1 or minus 1
        {
            //PointCollection textCoord_alternative = new PointCollection();
            if (CurrentMainOrdinalAngle == 179)
            {
                int z = 0;
            }
            
            if ((NumberOfDefreesToComplete_horizontal < 360) && (((CurrentMainOrdinalAngle + offset) * angleStep_horizontal) == NumberOfDefreesToComplete_horizontal))
           {
               return;
           }
            
            AngleVectorArray FirstCircle = VectorsArray[CurrentMainOrdinalAngle];
            AngleVectorArray SecondCircle = VectorsArray[(CurrentMainOrdinalAngle + offset) * angleStep_horizontal == NumberOfDefreesToComplete_horizontal ? 0 : (CurrentMainOrdinalAngle + offset)];
            bool flagOfend_horizontal = (CurrentMainOrdinalAngle + offset) * angleStep_horizontal == NumberOfDefreesToComplete_horizontal;
            AngleVector Firsti = new AngleVector();
            AngleVector Firstiplus1 = new AngleVector();
            AngleVector Secondi = new AngleVector();
            AngleVector Secondiplus1 = new AngleVector();
            bool continueflag;
            double marginfraction = 0;
            //double VerticalStepPerTextCoordinates = (1.0-marginfraction) / FirstCircle.array.Length;
            double VerticalStepPerTextCoordinates = (1.0 - marginfraction) / NumberOfDefreesToComplete_vertical / angleStep_vertical;
            double HorizontalStepPerTextCoordinates = (1.0-marginfraction) / (NumberOfDefreesToComplete_horizontal/angleStep_horizontal);
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

                
                if (i == FirstCircle.array.Length-1)
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
                    Firstiplus1 = (FirstCircle.array)[i+1];
                     Secondi=(SecondCircle.array)[i];
                     Secondiplus1 = (SecondCircle.array)[i+1];
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
                        double x = marginfraction/2.0+(double)CurrentMainOrdinalAngle *1.0* HorizontalStepPerTextCoordinates;
                        double y = marginfraction / 2.0 + (double)i * 1.0 * VerticalStepPerTextCoordinates;
                        
                         //  commented out
                        bool flagOfend_vertical= (Firsti.Angle == (360 - angleStep_vertical));
                        if (flagOfend_vertical == true)
                        {
                            int yyy = 0;
                        }
                        int lastcont = textCoord_alternative.Count;
                        textCoord_alternative.Add(new Point(x, y));
                        textCoord_alternative.Add(new Point(flagOfend_horizontal == true ? 1 : x + HorizontalStepPerTextCoordinates, y));
                        textCoord_alternative.Add(new Point(flagOfend_horizontal == true ? 1 : x + HorizontalStepPerTextCoordinates, flagOfend_vertical == true ? 1 : y + VerticalStepPerTextCoordinates));
                        textCoord_alternative.Add(new Point(x, flagOfend_vertical==true ? 1: y + VerticalStepPerTextCoordinates));
                       
                        
                        

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

        static public AngleVectorArray[] ReturnArrayOfAllPositionsInTube(int NumberOfDefreesToComplete_horizontal, int angleStep_horizontal, int NumberOfDefreesToComplete_vertical, int angleStep_vertical, double TubeThickness, double InternalDiameter, Vector3D AxisVector)
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
                Return360pointsGivenTwoPerPendicularVectors(NumberOfDefreesToComplete_vertical,angleStep_vertical, a_normalized, Final, ref arr, radiusPerPerpemdicularCircle, currFinalFinal);
                AngleVectorArray temp = new AngleVectorArray { Angle = i, Vector = currFinalFinal, array = arr };
                final[j] = temp;
                j++;
            }
            return final;

        }
        static private bool TryPerpendicularGivenVectorAndTwoPoints(Vector3D VectorInquestion, double x, double y, double z, ref Vector3D ResultVector)
        {

            bool result = true;
            double dotproductpartial;

            //if (PlaneInWhichPointIsSought == "Z")
            if (!(VectorInquestion.Z == 0))
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
            else if (!(VectorInquestion.Y == 0))// /(PlaneInWhichPointIsSought == "Y")
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
            if (result == true) ResultVector = new Vector3D(x, y, z);
            return result;

        }
        static private double FindAngleBtweemTwoVectors(Vector3D a, Vector3D b)
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
        static public void Return360pointsGivenTwoPerPendicularVectors(int NumberOfDefreesToComplete_vertical,int angleStep_vertical,Vector3D a, Vector3D ResultVector, ref AngleVector[] arr, double radiusPerPerpemdicularCircle, Vector3D CenterPointVector)
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
                AngleVector temp = new AngleVector { Angle = i, Vector = currFinalFinal };
                arr[j] = temp;
                j++;
                //dotProductTest = Vector3D.DotProduct(a, Final);
            }
        }
    }
}
