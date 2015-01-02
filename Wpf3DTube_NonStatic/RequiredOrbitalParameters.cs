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
    public class EarthOrbit
    {
        private Vector3D[] basis; // basis[0] is cosine leg and basis[1] is sine leg
        private Vector3D normalToEcliptic;
        private Vector3D vernalVector;
        private Vector3D periapsisdirection;
        private Vector3D tiltvector;
        public EarthOrbit()
        {
           
        }
        public Vector3D[] BASIS
        {
            get
            {
                return basis;
            }
            set
            {
                basis = value;
            }
        }
        public Vector3D EclipticNormal
        {
            get
            {
                return normalToEcliptic;
            }
            set
            {
                normalToEcliptic = value;
            }
        }
        public Vector3D VernalVector
        {
            get
            {
                return vernalVector;
            }
            set
            {
                vernalVector = value;
            }
        }
        public Vector3D PeriapsisDirection
        {
            get
            {
                return periapsisdirection;
            }
            set
            {
                periapsisdirection = value;
            }
        }

        public Vector3D TiltVector
        {
            get
            {
                return tiltvector;
            }
            set
            {
                tiltvector = value;
            }
        }

    }
    public class CosineandSineOfAngle
    {
        public double Sine { get; set; }
        public double Cosine { get; set; }
        public CosineandSineOfAngle()
        {

        }
        public CosineandSineOfAngle(double angle)
        {
            double PI = Math.PI;
            Sine = Math.Sin(angle / 180 * PI);
            Cosine = Math.Cos(angle / 180 * PI);
        }
    }
    public class RequiredOrbitalParameters
    {
        //public object epoch;
        public double aphelion;
        public double tilt;

        public double perihelion;

        public double semimajoraxis;

        public double eccentricity;
        public double inclination;

        public double longitudeoftheascendingnode;
        public double argumentofperiapsis;
        public object period;
        public string referencePlane;
        public RequiredOrbitalParameters()
        {
           
        }
        public RequiredOrbitalParameters(double asemimajoraxis, double aaphelion, double aperihelion, double atilt, double aeccentricity, double ainclination, double alongitudeoftheascendingnode, double aargumentofperiapsis, string areferencePlane)
        {
            aphelion = aaphelion;
            perihelion=aperihelion;
            tilt= atilt;
            eccentricity=aeccentricity;
            inclination=ainclination;
            longitudeoftheascendingnode= alongitudeoftheascendingnode;
            argumentofperiapsis=aargumentofperiapsis;
            referencePlane= areferencePlane;
            semimajoraxis = asemimajoraxis;
        }
        public EarthOrbit calculateOrbitRelatedVectors(PlanetPropeties sun, Vector3D normalOfPlanInWhichAscesnionNodeIS)
        {
            double Earthslongitudeoftheascendingnode = longitudeoftheascendingnode;
            Vector3D sunsAxisPerpendicularToItsEquator = sun.mormalToPlanetEquator;
            // ascension node direction is arbitrary in that it will determine the vernal vector and hence the point of periapsis. it is
            // arbitrary in virtue of providing arbitrary normal  (normalOfPlanInWhichAscesnionNodeIS) determining the plane in which it "resides"
     

            //GramSchmidt GM = new GramSchmidt(true, sunsAxisPerpendicularToItsEquator,h);
            Vector3D ascensiomModeVector = Vector3D.CrossProduct(sunsAxisPerpendicularToItsEquator,normalOfPlanInWhichAscesnionNodeIS);

            /// the normal to ecliptic MUST be perpendiculr to ascension node (bevasuse ascension node is on the ecliptic);ascension mde is also lies on sun's equator plane, so it's perpendicular to sun's axis.
            /// normal to eliptic is separated by 'inclination' number of degrees from that vector which is normal to sun's
            /// equator.  that means that cros product of ecliptic normal and sun's equator plane normal is the ascension node
            /// // such that with the provision that the angle between theose normal vectors is the inclination angle
            /// hence we use the folowing fucntion whic returmd the ecliptic normal (since we have one of the 2 normals, the angle between them and the vector which is perpendicular
            /// // to those normals)
            Vector3D eclipticNormal = ReturnVectorWithXAngleBetweenGivenVectorOnGivenPLane(inclination, sunsAxisPerpendicularToItsEquator, ascensiomModeVector);
            eclipticNormal.Normalize();

            Vector3D basisCosineLeg = ascensiomModeVector;
            Vector3D basisSineLeg = Vector3D.CrossProduct(eclipticNormal, ascensiomModeVector);
             basisSineLeg.Normalize();

            basisCosineLeg.Normalize();
            basisSineLeg.Normalize();
            CosineandSineOfAngle cosineandsine = new CosineandSineOfAngle(-360+Earthslongitudeoftheascendingnode);
            Vector3D vernalVector = basisCosineLeg * cosineandsine.Cosine + basisSineLeg * cosineandsine.Sine;
            vernalVector.Normalize();
            // direction of periapsis is argumentofperiaps is argumentofperiapsis degrees away from vernal vector (counter clock wise) on ecliptic
            Vector3D vectorOnEclipticPerpendicularToVernalVector = Vector3D.CrossProduct(eclipticNormal,vernalVector);
            vectorOnEclipticPerpendicularToVernalVector.Normalize();

            CosineandSineOfAngle cosineandsine_forargumentPeriapsis = new CosineandSineOfAngle(argumentofperiapsis);
            Vector3D vectorInDirectionInPeriapsis = vernalVector * cosineandsine_forargumentPeriapsis.Cosine + vectorOnEclipticPerpendicularToVernalVector * cosineandsine_forargumentPeriapsis.Sine;
            double focus = eccentricity * semimajoraxis;
            double distanceFromSunToEarthPeriapsis = semimajoraxis-focus;// focos distance
            Point3D periApsisPoint = sun.myCurrentLocation + vectorInDirectionInPeriapsis * distanceFromSunToEarthPeriapsis;
            Point3D imaginaryEarthOrbitCenter = periApsisPoint - semimajoraxis * vectorInDirectionInPeriapsis;
            Vector3D apoapsisDirectionVector=Vector3D.CrossProduct(eclipticNormal,vectorInDirectionInPeriapsis);
            apoapsisDirectionVector.Normalize();
            double semiminor = Math.Pow(Math.Pow(semimajoraxis, 2) - Math.Pow(focus, 2), 0.5);
            Point3D apoapsisPoint = sun.myCurrentLocation + apoapsisDirectionVector * semiminor;
            Vector3D tiltvector = ReturnVectorWithXAngleBetweenGivenVectorOnGivenPLane(tilt, eclipticNormal, vernalVector);
            EarthOrbit EO = new EarthOrbit();
            EO.BASIS = new Vector3D[] { vectorInDirectionInPeriapsis, apoapsisDirectionVector };
            EO.EclipticNormal = eclipticNormal;
            EO.VernalVector = vernalVector;
            EO.TiltVector = tiltvector;
            EO.PeriapsisDirection = vectorInDirectionInPeriapsis;
            return EO;

        }
        public static Vector3D calculateVernalPoint(PlanetPropeties earthPlanet, Vector3D  normalOfPlaneInWhichTiltVectorProducesTiltAngleWithYAxix)
        {
            // normalOfPlaneInWhichTiltVectorProducesTiltAngleWithYAxix is normally   new Vector3D(0, 1, 0);
            double tilt = earthPlanet.Tilt;
            Vector3D axisIfTiltisZero = new Vector3D(0, 1, 0);
            
            Vector3D ExlipticPlanev1 = normalOfPlaneInWhichTiltVectorProducesTiltAngleWithYAxix;
            Vector3D ExlipticPlanev2 = Vector3D.CrossProduct(axisIfTiltisZero, ExlipticPlanev1);
            Vector3D vectorRepresentingRealAxis = ReturnVectorWithXAngleBetweenGivenVectorOnGivenPLane(tilt, axisIfTiltisZero, normalOfPlaneInWhichTiltVectorProducesTiltAngleWithYAxix);

            /// we aaume that the north pole points in the direction of y axis
            /// // the equatorial plane is the xz palne
            /// the tilt is 23.5 degrees projecting on x axis
            Vector3D resultVector;
            resultVector = Vector3D.CrossProduct(axisIfTiltisZero, vectorRepresentingRealAxis);
            resultVector.Normalize();
            return resultVector;
        }
        public static Vector3D ReturnVectorOffPeriapsis(double argumanetOfPeriapsis,Vector3D ascensionNodeVector,Vector3D normalToOrvitPlane)
        {
            Vector3D perpendicularToAscensionNodeAndOnPlane = Vector3D.CrossProduct(normalToOrvitPlane, ascensionNodeVector);
            perpendicularToAscensionNodeAndOnPlane.Normalize();
            /// now use perpendicularToAscensionNodeAndOnPlane and ascensionNodeVector as basis n order to rotate ascensop node counter clockwise by argumanetOfPeriapsis degrees



            double PI = Math.PI;
            double sinA = Math.Sin(argumanetOfPeriapsis / 180 * PI);
            double cosA = Math.Cos(argumanetOfPeriapsis / 180 * PI);
            Vector3D result = cosA * ascensionNodeVector + sinA * perpendicularToAscensionNodeAndOnPlane;
            return result;
        }

        static private Vector3D ReturnVectorWithXAngleBetweenSaidVectorAndPlane(double x, Vector3D vectorToWhichInclinationRefersTo, Vector3D normal)
        {
            //Vector3D vectorPerPendicularToBothGivenVectors = Vector3D.CrossProduct(vectorToWhichInclinationRefersTo, normal);
            //vectorPerPendicularToBothGivenVectors = vectorPerPendicularToBothGivenVectors/vectorPerPendicularToBothGivenVectors.Length;
            //// we have the sun's  normal to its equator plane plus a vector which is on that plane (preoendicular to the normal) and which is a projection
            /// of the major axis vecotr with magnitude of 1
            // now we find the vector which will be on the plane of the body in question taking into axccount the inclination
            double PI = Math.PI;
            double sinA = Math.Sin(x / 180 * PI);
            double cosA = Math.Cos(x / 180 * PI);
            // weassume the normal and the projection vector are unti vectors;
            Vector3D v1 = normal * sinA;
            Vector3D v2 = vectorToWhichInclinationRefersTo * cosA;
            return (v1 + v2) / (v1 + v2).Length;

        }

        static public Vector3D ReturnVectorWithXAngleBetweenGivenVectorOnGivenPLane(double x, Vector3D firstLegOfTheAngle, Vector3D normalVorerspondingToDesierdPlane)
        {

            Vector3D n = normalVorerspondingToDesierdPlane;
            n.Normalize();
            Vector3D v1 = firstLegOfTheAngle;
            // v2 is perpendiculatr to firstLegOfTheAngle and is pn the same plane relating to our given normal vector; thus it is the second leg of the basis for the plane which
            // is perpendiculat to the vector
            Vector3D v2 = Vector3D.CrossProduct(v1, n); ;
            v2.Normalize();
            Vector3D resultVector;
            resultVector = RotateBasisByXDegrees(x, v1, v2); ;
            return resultVector;
        }
        public static Vector3D RotateBasisByXDegrees(double x,Vector3D v1_cosdirection,Vector3D v2_sindirection)
        {
            Vector3D resultVector;
            double PI = Math.PI;
            double sinA = Math.Sin(x / 180 * PI);
            double cosA = Math.Cos(x / 180 * PI);
            resultVector = v1_cosdirection * cosA + v2_sindirection * sinA;
            return resultVector;
        }

    }
}
