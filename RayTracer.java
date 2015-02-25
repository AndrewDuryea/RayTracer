/*
 * Andrew Duryea
 * February 25, 2013
 * RayTracer.java
 *
 * This program is a simple ray tracer.
 */

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.border.TitledBorder;
import java.util.ArrayList;

public class RayTracer extends JFrame
{
	//JPanel needs to be subclasses to override the paintComponent method
	JPanel canvas = new JPanel()
	{
		public void paintComponent(Graphics g)
		{
			super.paintComponent(g);
			//this will also move the canvas, but is too easy
			//g.translate(canvas.getWidth() / 2,canvas.getHeight() / 2);

			double xOffset = canvas.getWidth() / 2;
			double yOffset = canvas.getHeight() / 2;

			Vector3D intensity;
			Vector3D ray;

			for(double x = -xOffset; x < xOffset; x++)
			{
				for(double y = -yOffset; y < yOffset; y++)
				{
					ray = new Vector3D(x - centerOfProjection.X, y - centerOfProjection.Y, 0 - centerOfProjection.Z);

					intensity = traceRay(MAXDEPTH, centerOfProjection, ray, 1, true);

					intensity.end.X = (intensity.end.X > 1) ? 1 : intensity.end.X;
					intensity.end.Y = (intensity.end.Y > 1) ? 1 : intensity.end.Y;
					intensity.end.Z = (intensity.end.Z > 1) ? 1 : intensity.end.Z;

					intensity.end.X = (intensity.end.X < 0) ? 1 : intensity.end.X;
					intensity.end.Y = (intensity.end.Y < 0) ? 1 : intensity.end.Y;
					intensity.end.Z = (intensity.end.Z < 0) ? 1 : intensity.end.Z;

					g.setColor(new Color((int)(255*intensity.end.X),(int)(255*intensity.end.Y),(int)(255*intensity.end.Z)));
					g.drawLine((int)(x+xOffset),(int)(-y+yOffset),(int)(x+xOffset),(int)(-y+yOffset));
				}
			}
		}
	};

	//the selected light
	PointLight currentLight;

	double Iar = 0.5, Iag = 0.5, Iab = 0.5; //the ambient intensity
	ArrayList<PointLight> pointLightList = new ArrayList<PointLight>();

	ArrayList<Object3D> objectList = new ArrayList<Object3D>();

	static int windowHeight = 650;
	static int windowWidth = 850;

	Point3D centerOfProjection = new Point3D(0,0,-800);

	double threshold = 0.01; //threshold for adaptive depth

	double airIndex = 1;
	int MAXDEPTH = 5;
	double intersectDist;
	Object3D currentObject;
	Point3D intersection;
	Vector3D normal;

	Vector3D back = new Vector3D(0.6,0.6,1); //the background intensity

	public RayTracer()
	{
		objectList.add(new Sphere(-200,-100,-100,100,1,0,0,0.2,0,0.8));	//sphere 1
		objectList.add(new Sphere(-100,0,500,200,0,.5,0,0.5,0.5,0));	//sphere 2
		objectList.add(new Sphere(310,0,500,200,0,0,.5,0.5,0.5,0));	//sphere 3
		objectList.add(new Sphere(0,-450,500,200,0.78,0.5,0,0.5,0,0));	//sphere 4
		objectList.add(new Checkerboard(0,-250,0,100,100));//checkerboard

		currentLight = new PointLight(70,140,140,140, new Vector3D(120,210,-150));
		pointLightList.add(currentLight);

		this.setLayout(new BorderLayout());
			canvas.setLayout(new FlowLayout());

		add(canvas, BorderLayout.CENTER);

		canvas.setFocusable(true);
		canvas.requestFocusInWindow();

		//canvas.setBackground(Color.white);
		canvas.setBackground(new Color(0,0,0));

    }

    public static void main(String s[])
    {
    	RayTracer f = new RayTracer();

    	//window listener so close works
		f.addWindowListener(new WindowAdapter()
        {
            public void windowClosing(WindowEvent e) {System.exit(0);}
        });

		f.setBounds(200, 50, windowWidth, windowHeight);
		f.setTitle("RayTracer");
		f.setVisible(true);
    }

    public Vector3D traceRay(int depth, Point3D projectP, Vector3D ray, double reflecConst, boolean isEntering)
    {
    	//System.out.println("here " + depth);
		if(depth <= 0)
		{
			//System.out.println("here");
			return new Vector3D(0,0,0);
		}

		intersectDist = 100000;
		currentObject = null;

		for(Object3D o : objectList)
		{
			o.getIntersection(projectP, ray);
		}

		if(currentObject != null)
		{
			return currentObject.getIntensity(depth, intersection, ray, reflecConst, isEntering);
		}
		else
		{
    		return back; //the background color
		}
    }


    public boolean shadowTraceRay(Point3D projectP, Vector3D ray, Vector3D shadow, Double dist)
    {
		intersectDist = 100000;
		currentObject = null;

		for(Object3D o : objectList)
		{
			o.getIntersection(projectP, ray);
		}

		if(currentObject != null)
		{
			Vector3D s = currentObject.isTranslucent();

    		if(s != null)
	    		shadow.end = s.end;

	    	dist = intersectDist;

			return true;
		}
		else
		{
    		return false;
		}
    }

	boolean project = true;

	//the default veiwing distance
	double viewpoint = 200;

	//the view vector
	Vector3D view = new Vector3D(0,0,1);

	ArrayList<Object3D> objects = new ArrayList<Object3D>();

	//this class is a simple point in 3D
	public class Point3D
	{
		double X;
		double Y;
		double Z;

		String pointID;

		Point3D(double X, double Y, double Z)
		{
			this.X = X;
			this.Y = Y;
			this.Z = Z;
			this.pointID = ""+-1; //this point is unlisted, ie not part of a triangle or temporary
		}

		Point3D(double X, double Y, double Z, int pointID)
		{
			this.X = X;
			this.Y = Y;
			this.Z = Z;
			this.pointID = "P" + pointID; //this point is listed
		}

		Point3D(double X, double Y, double Z, String name)
		{
			this.X = X;
			this.Y = Y;
			this.Z = Z;
			this.pointID = name; //this is a named point
		}

		public void setPoint(double X, double Y, double Z)
		{
			this.X = X;
			this.Y = Y;
			this.Z = Z;
		}

		Point3D perspectiveProject(double view)
		{
			double Xps;
			double Yps;
			double Zps;

			if(project)
			{
				Xps = (view * X) / (view + Z);
				Yps = (view * Y) / (view + Z);
				Zps = (Z) / (view + Z);
			}
			else
			{
				Xps = X;
				Yps = Y;
				Zps = Z;
			}
			//System.out.println("here3" + faceList.size());
			return new Point3D(Xps, Yps, Zps);
		}

/*		public Point3D multiply(Matrix m)
		{
			Matrix point = new Matrix(1,4);

			//create a matrix from the point
			point.setValue(0,0,this.X);
			point.setValue(0,1,this.Y);
			point.setValue(0,2,this.Z);
			point.setValue(0,3,1);

			//pull the point through the composite matrix
			point = point.multiply(m);

			return new Point3D(point.getValue(0,0), point.getValue(0,1), point.getValue(0,2));
		}
*/
		public String toString()
		{
			return "(" + this.X + ", " + this.Y + ", " + this.Z + ")";
		}

		//lightStuff
		/*public Vector3D getIntensityVector(Vector3D normal, Object3D o, Point3D point)
		{
			//the reflection vector
    		Vector3D Ref = null;

    		double cphi;
    		double ctheta;

    		double dist;

			double Ir, Ig, Ib;
			Ir = Ig = Ib = 0;

			//ambient diffuse
			Ir = Iar * o.Kdr;
			Ig = Iag * o.Kdg;
			Ib = Iab * o.Kdb;

			for(PointLight p : pointLightList)
			{
				cphi = normal.unitVector().dot(p.direction.unitVector());
				dist = p.distance;

				//specular
				if(cphi > 0) //case 1
				{

					Ir += p.Ipr * o.Kdr * cphi / dist;
					Ig += p.Ipg * o.Kdg * cphi / dist;
					Ib += p.Ipb * o.Kdb * cphi / dist;

					Ref = normal.unitVector().minus(p.direction.unitVector().divideScaler(2*cphi));

					Ir += p.Ipr / dist * o.Ks * Math.pow((Ref.unitVector().dot(currentCamera.viewVec.unitVector())), o.n);
					Ig += p.Ipg / dist * o.Ks * Math.pow((Ref.unitVector().dot(currentCamera.viewVec.unitVector())), o.n);
					Ib += p.Ipb / dist * o.Ks * Math.pow((Ref.unitVector().dot(currentCamera.viewVec.unitVector())), o.n);
				}
				else //case 3
				{
					Ref = p.direction.unitVector().flipDirection();

					Ir += p.Ipr / dist * o.Ks * Math.pow((Ref.unitVector().dot(currentCamera.viewVec.unitVector())), o.n);
					Ig += p.Ipg / dist * o.Ks * Math.pow((Ref.unitVector().dot(currentCamera.viewVec.unitVector())), o.n);
					Ib += p.Ipb / dist * o.Ks * Math.pow((Ref.unitVector().dot(currentCamera.viewVec.unitVector())), o.n);
				}
			}

			Ir = (Ir > 1) ? 1 : Ir;
			Ig = (Ig > 1) ? 1 : Ig;
			Ib = (Ib > 1) ? 1 : Ib;

			Ir = (Ir < 0) ? 0 : Ir;
			Ig = (Ig < 0) ? 0 : Ig;
			Ib = (Ib < 0) ? 0 : Ib;

			return new Vector3D(Ir, Ig, Ib);
		}*/
	}

    //this class is a vector in 3D
	public class Vector3D
	{
		Point3D start = new Point3D(0,0,0);
		Point3D end;

		//find the vector of two points
		Vector3D(Point3D start, Point3D end)
		{
			this.end = new Point3D(end.X - start.X, end.Y-start.Y, end.Z-start.Z);
		}

		//for when we know where the vector ends
		Vector3D(double X, double Y, double Z)
		{
			this.end = new Point3D(X,Y,Z);
		}

		//for when we know where the vector ends
		Vector3D(Point3D end)
		{
			this.end = end;
		}

		public Vector3D add(Vector3D Q)
		{
			//System.out.println(this.end.X + " " + Q.end.X);
			return new Vector3D(this.end.X+Q.end.X, this.end.Y+Q.end.Y, this.end.Z+Q.end.Z);
		}

		public Vector3D minus(Vector3D Q)
		{
			return new Vector3D(this.end.X-Q.end.X, this.end.Y-Q.end.Y, this.end.Z-Q.end.Z);
		}

		public Vector3D addScaler(double c)
		{
			return new Vector3D(this.end.X+c, this.end.Y+c, this.end.Z+c);
		}

		public Vector3D minusScaler(double c)
		{
			return new Vector3D(this.end.X-c, this.end.Y-c, this.end.Z-c);
		}

		public Vector3D divideScaler(double c)
		{
			return new Vector3D(this.end.X/c, this.end.Y/c, this.end.Z/c);
		}

		public Vector3D multiplyScaler(double c)
		{
			return new Vector3D(this.end.X*c, this.end.Y*c, this.end.Z*c);
		}

		public Vector3D flipDirection()
		{
			return new Vector3D(-this.end.X,-this.end.Y,-this.end.Z);
		}

		//the cross product of two vectors, this*Q
		public Vector3D cross(Vector3D Q)
		{
			//this	x y z
			//Q		x y z
			return new Vector3D(
				this.end.Y*Q.end.Z-this.end.Z*Q.end.Y,
				-(this.end.X*Q.end.Z-this.end.Z*Q.end.X),
				this.end.X*Q.end.Y-this.end.Y*Q.end.X
			);
		}

		//the dot product of two vectors, this.Q
		public double dot(Vector3D Q)
		{
			return this.end.X*Q.end.X+this.end.Y*Q.end.Y+this.end.Z*Q.end.Z;
		}

		public double magnitude()
		{
			return Math.sqrt(end.X*end.X + end.Y*end.Y + end.Z*end.Z);
		}

		public Vector3D unitVector()
		{
			return new Vector3D(end.X/this.magnitude(),end.Y/this.magnitude(),end.Z/this.magnitude());
		}

		public String toString()
		{
			return "<" + end.X + ", " + end.Y + ", " + end.Z + ">";
		}

		public Vector3D perspectiveProject(double viewpoint)
		{
			return new Vector3D(end.perspectiveProject(viewpoint));
		}
	}

	public interface Object3D
	{
		public void getIntersection(Point3D projectP, Vector3D ray);
		public Vector3D getIntensity(int depth, Point3D projectP, Vector3D ray, double reflecConst, boolean isEntering);
		public Vector3D isTranslucent();
		public double getTransmitted();
	}

	public class Sphere implements Object3D
	{
		double X;
		double Y;
		double Z;

		double radius;

		double Kdr, Kdg, Kdb;

		double Ks = .5;
		double n = 40;

		double local, reflected, transmitted;

		double refracIndex = 1.5; //~index of glass

		public Sphere(double X, double Y, double Z, double radius, double Kdr, double Kdg, double Kdb, double local, double reflected, double transmitted)
		{
			this.X = X;
			this.Y = Y;
			this.Z = Z;
			this.radius = radius;
			this.Kdr = Kdr;
			this.Kdg = Kdg;
			this.Kdb = Kdb;
			this.local = local;
			this.reflected = reflected;
			this.transmitted = transmitted;
		}

		public void getIntersection(Point3D projectP, Vector3D ray)
		{
			double h = (projectP.X-this.X)*(projectP.X-this.X) + (projectP.Y-this.Y)*(projectP.Y-this.Y) + (projectP.Z-this.Z)*(projectP.Z-this.Z) - radius*radius;
			if(Math.abs(h) <= 0.001)
				return; //if the point is on the sphere it will always intersect

			double disc, ts1, ts2;
			double l, m, n, r, asphere, bsphere, csphere, tsphere;

			asphere = ray.end.X *ray.end.X + ray.end.Y * ray.end.Y + ray.end.Z * ray.end.Z;
			bsphere = 2 * ray.end.X * (projectP.X - this.X) + 2 * ray.end.Y * (projectP.Y - this.Y) + 2 * ray.end.Z * (projectP.Z - this.Z);
			csphere = this.X * this.X + this.Y * this.Y + this.Z * this.Z + projectP.X * projectP.X + projectP.Y * projectP.Y + projectP.Z * projectP.Z +
				2 * (-this.X * projectP.X - this.Y * projectP.Y - this.Z * projectP.Z) - radius * radius;

			disc = bsphere * bsphere - 4 * asphere * csphere;

			if (disc < 0)
				return; //no intersection

			ts1 = (-bsphere + Math.sqrt(disc)) / (2 * asphere);
			ts2 = (-bsphere - Math.sqrt(disc)) / (2 * asphere);
			if(ts1 >= ts2)
				tsphere = ts2;
			else
				tsphere = ts1;

			if(intersectDist < tsphere)
				return; //another object is closer

			if(tsphere < 0)
				return; //no visible intersection

			intersectDist = tsphere;
			currentObject = this;
			intersection = new Point3D((projectP.X + ray.end.X * tsphere),(projectP.Y + ray.end.Y * tsphere),(projectP.Z + ray.end.Z * tsphere));
			normal = new Vector3D((intersection.X-this.X), (intersection.Y-this.Y), (intersection.Z-this.Z));

			return;
		}

		public Vector3D getIntensity(int depth, Point3D projectP, Vector3D ray, double reflecConst, boolean isEntering)
		{
			Vector3D ret = new Vector3D(0,0,0);

			Vector3D shadow = new Vector3D(0,0,0);

			Double shadowDist;

			for(PointLight p : pointLightList)
			{
				shadowDist = p.distance;

				if(shadowTraceRay(projectP, p.direction, shadow, shadowDist))
				{
					ret = new Vector3D(local*Kdr, local*Kdg, local*Kdb);

					return ret.add(shadow.multiplyScaler(shadowDist.doubleValue()/p.distance));
				}
			}

			//if((reflecConst * reflected <= threshold)||(reflecConst * transmitted <= threshold))
			//	return back;

			Vector3D reflec;
			Vector3D intensity;
			Vector3D normalUnit = normal.unitVector();
			Vector3D rayUnit = ray.unitVector();

			if(local >= 0.0)
			{
				double cphi;
				double ctheta;

				Vector3D Ref;

	    		double dist;

				double Ir, Ig, Ib;
				Ir = Ig = Ib = 0;

				//ambient diffuse
				Ir = Iar * Kdr;
				Ig = Iag * Kdg;
				Ib = Iab * Kdb;

				for(PointLight p : pointLightList)
				{
					cphi = normal.unitVector().dot(p.direction.unitVector());
					dist = p.distance;
					//dist = Math.sqrt(Math.pow((p.direction.end.X - intersection.X),2) + Math.pow((p.direction.end.Y - intersection.Y),2) + Math.pow((p.direction.end.Z - intersection.Z),2)) + p.distance;

					if(cphi > 0) //case 1
					{
						//System.out.println("here1");

						//diffuse
						Ir += p.Ipr * Kdr * cphi / dist;
						Ig += p.Ipg * Kdg * cphi / dist;
						Ib += p.Ipb * Kdb * cphi / dist;

						Ref = normal.unitVector().minus(p.direction.unitVector().divideScaler(2*cphi));

						//specular
						Ir += p.Ipr / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
						Ig += p.Ipg / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
						Ib += p.Ipb / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
					}
					else //case 3
					{
						//System.out.println("here2");
						Ref = p.direction.unitVector().flipDirection();

						//specular
						Ir += p.Ipr / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
						Ig += p.Ipg / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
						Ib += p.Ipb / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
					}
				}

				ret = ret.add(new Vector3D((local * Ir),(local * Ig),(local * Ib)));
			}

			if(reflected >= 0.0)
			{
				double cphi = (-rayUnit.end.X) * normalUnit.end.X +
					(-rayUnit.end.Y) * normalUnit.end.Y +
					(-rayUnit.end.Z) * normalUnit.end.Z;

				if(cphi > 0)
				{
					reflec = new Vector3D((normalUnit.end.X - (-rayUnit.end.X) / (2 * cphi)),
						(normalUnit.end.Y - (-rayUnit.end.Y) / (2 * cphi)),
						(normalUnit.end.Z - (-rayUnit.end.Z) / (2 * cphi)));
				}
				else if(cphi == 0)
				{
					reflec = rayUnit;
				}
				else
				{
					reflec = new Vector3D((-normalUnit.end.X + (-rayUnit.end.X) / (2 * cphi)),
						(-normalUnit.end.Y + (-rayUnit.end.Y) / (2 * cphi)),
						(-normalUnit.end.Z + (-rayUnit.end.Z) / (2 * cphi)));
				}

				intensity = traceRay((depth-1), projectP, reflec, (reflecConst * reflected), true);

				intensity = intensity.multiplyScaler(reflected);
				//System.out.println(intensity);

				ret = ret.add(intensity);
			}

			if(transmitted >= 0.0)
			{
				double nu;

				//if(isEntering)
				//{
					nu = refracIndex / airIndex;
				//}
				//else
				//{
				//	nu = airIndex / refracIndex;
				//}

				Vector3D a = ray.unitVector().divideScaler(nu);
				double b = Math.sqrt(1-(1/(nu*nu))*(1-Math.pow((ray.unitVector().flipDirection().dot(normal.unitVector())),2)));
				//double b = Math.sqrt(1-(1/(nu*nu))*(1-(ray.unitVector().flipDirection().dot(normal.unitVector()))));
				double c = (1/nu) * (ray.unitVector().flipDirection().dot(normal.unitVector()));

				Vector3D refrac = a.minus(normal.unitVector().multiplyScaler(b-c));

				ret = ret.add(traceRay((depth-1), projectP, refrac, (reflecConst * transmitted), !isEntering).multiplyScaler(transmitted));

				//ret = ret.add(traceRay((depth-1), projectP, ray, (reflecConst * transmitted)).multiplyScaler(transmitted));
			}

			return ret;

		}

		public Vector3D isTranslucent()
		{
			if(transmitted >= 0)
			{
				return new Vector3D(transmitted * Kdr, transmitted * Kdg, transmitted * Kdb);
			}

			return null;
		}

		public double getTransmitted()
		{
			return transmitted;
		}
	}

	public class Checkerboard implements Object3D
	{
		double X;
		double Y;
		double Z;

		double height;
		double width;

		double Ks = .5;
		double n = 10;

		double local = 0.6;
		double reflected = 0.4;
		double reflected2 = 0.1;
		double transmitted = 0.3;

		double refracIndex = 1.2;

		public Checkerboard(double x, double y, double z, double height, double width)
		{
			this.X = x;
			this.Y = y;
			this.Z = z;
			this.height = height;
			this.width = width;
		}

		public void getIntersection(Point3D projectP, Vector3D ray)
		{
			double a, b, c;
			double denom, objectDist, d, x, y, z;

			//normal <a,b,c> of plane
			a = 0;
			b = 1;
			c = 0;

			Vector3D n = new Vector3D(a,b,c);

			Vector3D P1 = new Vector3D(X,Y,Z);
			Vector3D P2 = new Vector3D(projectP);

			if(Math.abs(n.dot(P2.minus(P1))) <= 0.001)
				return; //the point is in the plane, and will thus always intersect

			denom = a * ray.end.X + b * ray.end.Y + c * ray.end.Z;

			if(Math.abs(denom) <= 0.001)
				return; //ray parallel to plane

			d = a * this.X + b * this.Y + c * this.Z;

			objectDist = -(a * projectP.X + b * projectP.Y + c * projectP.Z - d) / denom;
			x = projectP.X + ray.end.X * objectDist;
			y = projectP.Y + ray.end.Y * objectDist;
			z = projectP.Z + ray.end.Z * objectDist;

			if((z < 0) || (z > 1000) || (objectDist < 0)) //no visible intersection
			{
				return;
			}

			if((x < -500) || (x > 500) || (objectDist < 0)) //no visible intersection
			{
				return;
			}

			if(intersectDist < objectDist) //another object is nearer
			{
				return;
			}

			//System.out.println(ray);
			//System.out.println();

			intersectDist = objectDist;
			currentObject = this;
			intersection = new Point3D(x, y, z);
			normal = new Vector3D(0,1,0);
		}

		public Vector3D getIntensity(int depth, Point3D projectP, Vector3D ray, double reflecConst, boolean isEntering)
		{
			boolean bColor;

			//Vector3D normal = new Vector3D(0,1,0);

			double Kdr = 0, Kdg = 0, Kdb = 0;

			double Ir, Ig, Ib;
			Ir = Ig = Ib = 0;

			double cphi;
			double ctheta;

			Vector3D Ref;

    		double dist;

			Vector3D reflec;
			Vector3D intensity;
			Vector3D normalUnit = normal.unitVector();
			Vector3D rayUnit = ray.unitVector();
			Vector3D ret = new Vector3D(0,0,0);
			Vector3D shadow  = new Vector3D(0,0,0);

			Double shadowDist;

			if(intersection.X >= 0)
				bColor = true;
			else
				bColor = false;

			if(Math.abs(intersection.X%400) > 200)
				bColor = !bColor;
			if(Math.abs(intersection.Z%400) > 200)
				bColor = !bColor;

			if(bColor) //red
			{
				Kdr = 1;
				Kdg = 0;
				Kdb = 0;


				for(PointLight p : pointLightList)
				{
					shadowDist = p.distance;

					if(shadowTraceRay(projectP, p.direction, shadow, shadowDist))
					{
						ret = new Vector3D(local*Kdr, local*Kdg, local*Kdb);

						return ret.add(shadow.multiplyScaler(shadowDist.doubleValue()/p.distance));
					}
				}

				Ir = Ig = Ib = 0;

				if(local >= 0.0)
				{
					//ambient diffuse
					Ir = Iar * Kdr;
					Ig = Iag * Kdg;
					Ib = Iab * Kdb;

					for(PointLight p : pointLightList)
					{
						cphi = normal.unitVector().dot(p.direction.unitVector());
						dist = p.distance;
						//dist = Math.sqrt(Math.pow((p.direction.end.X - intersection.X),2) + Math.pow((p.direction.end.Y - intersection.Y),2) + Math.pow((p.direction.end.Z - intersection.Z),2)) + p.distance;

						if(cphi > 0) //case 1
						{
							//System.out.println("here1");

							//diffuse
							Ir += p.Ipr * Kdr * cphi / dist;
							Ig += p.Ipg * Kdg * cphi / dist;
							Ib += p.Ipb * Kdb * cphi / dist;

							Ref = normal.unitVector().minus(p.direction.unitVector().divideScaler(2*cphi));

							//specular
							Ir += p.Ipr / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
							Ig += p.Ipg / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
							Ib += p.Ipb / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
						}
						else //case 3
						{
							//System.out.println("here2");
							Ref = p.direction.unitVector().flipDirection();

							//specular
							Ir += p.Ipr / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
							Ig += p.Ipg / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
							Ib += p.Ipb / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
						}
					}

					ret = ret.add(new Vector3D((local * Ir),(local * Ig),(local * Ib)));
				}

				if(reflected >= 0.0)
				{
					cphi = (-rayUnit.end.X) * normalUnit.end.X +
						(-rayUnit.end.Y) * normalUnit.end.Y +
						(-rayUnit.end.Z) * normalUnit.end.Z;

					if(cphi > 0)
					{
						reflec = new Vector3D((normalUnit.end.X - (-rayUnit.end.X) / (2 * cphi)),
							(normalUnit.end.Y - (-rayUnit.end.Y) / (2 * cphi)),
							(normalUnit.end.Z - (-rayUnit.end.Z) / (2 * cphi)));
					}
					else if(cphi == 0)
					{
						reflec = rayUnit;
					}
					else
					{
						reflec = new Vector3D((-normalUnit.end.X + (-rayUnit.end.X) / (2 * cphi)),
							(-normalUnit.end.Y + (-rayUnit.end.Y) / (2 * cphi)),
							(-normalUnit.end.Z + (-rayUnit.end.Z) / (2 * cphi)));
					}

					intensity = traceRay((depth-1), projectP, reflec, (reflecConst * reflected), true);

					intensity = intensity.multiplyScaler(reflected);
					//System.out.println(intensity);

					ret = ret.add(intensity);
				}
			}
			else //white
			{
				Kdr = 0.5;
				Kdg = 0.5;
				Kdb = 0.5;


				for(PointLight p : pointLightList)
				{
					shadowDist = p.distance;

					if(shadowTraceRay(projectP, p.direction, shadow, shadowDist))
					{
						ret = new Vector3D(local*Kdr, local*Kdg, local*Kdb);

						return ret.add(shadow.multiplyScaler(shadowDist.doubleValue()/p.distance));
					}
				}

				Ir = Ig = Ib = 0;

				if(local >= 0.0)
				{
					//ambient diffuse
					Ir = Iar * Kdr;
					Ig = Iag * Kdg;
					Ib = Iab * Kdb;

					for(PointLight p : pointLightList)
					{
						cphi = normal.unitVector().dot(p.direction.unitVector());
						dist = p.distance;
						//dist = Math.sqrt(Math.pow((p.direction.end.X - intersection.X),2) + Math.pow((p.direction.end.Y - intersection.Y),2) + Math.pow((p.direction.end.Z - intersection.Z),2)) + p.distance;

						if(cphi > 0) //case 1
						{
							//System.out.println("here1");

							//diffuse
							Ir += p.Ipr * Kdr * cphi / dist;
							Ig += p.Ipg * Kdg * cphi / dist;
							Ib += p.Ipb * Kdb * cphi / dist;

							Ref = normal.unitVector().minus(p.direction.unitVector().divideScaler(2*cphi));

							//specular
							Ir += p.Ipr / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
							Ig += p.Ipg / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
							Ib += p.Ipb / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
						}
						else //case 3
						{
							//System.out.println("here2");
							Ref = p.direction.unitVector().flipDirection();

							//specular
							Ir += p.Ipr / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
							Ig += p.Ipg / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
							Ib += p.Ipb / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
						}
					}

					ret = ret.add(new Vector3D((local * Ir),(local * Ig),(local * Ib)));
				}

				if(reflected2 >= 0.0)
				{
					cphi = (-rayUnit.end.X) * normalUnit.end.X +
						(-rayUnit.end.Y) * normalUnit.end.Y +
						(-rayUnit.end.Z) * normalUnit.end.Z;

					if(cphi > 0)
					{
						reflec = new Vector3D((normalUnit.end.X - (-rayUnit.end.X) / (2 * cphi)),
							(normalUnit.end.Y - (-rayUnit.end.Y) / (2 * cphi)),
							(normalUnit.end.Z - (-rayUnit.end.Z) / (2 * cphi)));
					}
					else if(cphi == 0)
					{
						reflec = rayUnit;
					}
					else
					{
						reflec = new Vector3D((-normalUnit.end.X + (-rayUnit.end.X) / (2 * cphi)),
							(-normalUnit.end.Y + (-rayUnit.end.Y) / (2 * cphi)),
							(-normalUnit.end.Z + (-rayUnit.end.Z) / (2 * cphi)));
					}

					intensity = traceRay((depth-1), projectP, reflec, (reflecConst * reflected2), true);

					intensity = intensity.multiplyScaler(reflected2);
					//System.out.println(intensity);

					ret = ret.add(intensity);
				}

				if(transmitted >= 0.0)
				{
					double nu;

					//if(isEntering)
					//{
						nu = refracIndex / airIndex;
					//}
					//else
					//{
					//	nu = airIndex / refracIndex;
					//}

					Vector3D a = ray.unitVector().divideScaler(nu);
					double b = Math.sqrt(1-(1/(nu*nu))*(1-Math.pow((ray.unitVector().flipDirection().dot(normal.unitVector())),2)));
					//double b = Math.sqrt(1-(1/(nu*nu))*(1-(ray.unitVector().flipDirection().dot(normal.unitVector()))));
					double c = (1/nu) * (ray.unitVector().flipDirection().dot(normal.unitVector()));

					Vector3D refrac = a.minus(normal.unitVector().multiplyScaler(b-c));

					ret = ret.add(traceRay((depth-1), projectP, refrac, (reflecConst * transmitted), !isEntering).multiplyScaler(transmitted));

					//ret = ret.add(traceRay((depth-1), projectP, ray, (reflecConst * transmitted), !isEntering).multiplyScaler(transmitted));
				}
			}


			return ret;

		}

/*
		public Vector3D getIntensity(int depth, Point3D projectP, Vector3D ray, double reflecConst)
		{
			boolean bColor;

			Vector3D normal = new Vector3D(0,0,1);

			double Kdr = 0, Kdg = 0, Kdb = 0;

			if(intersection.X >= 0)
				bColor = true;
			else
				bColor = false;

			if(Math.abs(intersection.X%400) > 200)
				bColor = !bColor;
			if(Math.abs(intersection.Z%400) > 200)
				bColor = !bColor;

			if(bColor) //red
			{
				Kdr = 1;
				Kdg = 0;
				Kdb = 0;
			}
			else //white
			{
				Kdr = 1;
				Kdg = 1;
				Kdb = 1;
			}

			//the reflection vector
    		Vector3D Ref = null;

    		double cphi;
    		double ctheta;

    		double dist;

			double Ir, Ig, Ib;
			Ir = Ig = Ib = 0;

			//ambient diffuse
			Ir = Iar * Kdr;
			Ig = Iag * Kdg;
			Ib = Iab * Kdb;

			for(PointLight p : pointLightList)
			{
				cphi = normal.unitVector().dot(p.direction.unitVector());
				dist = p.distance;
				//dist = Math.sqrt(Math.pow((p.direction.end.X - intersection.X),2) + Math.pow((p.direction.end.Y - intersection.Y),2) + Math.pow((p.direction.end.Z - intersection.Z),2)) + p.distance;

				if(cphi > 0) //case 1
				{
					//System.out.println("here1");

					//diffuse
					Ir += p.Ipr * Kdr * cphi / dist;
					Ig += p.Ipg * Kdg * cphi / dist;
					Ib += p.Ipb * Kdb * cphi / dist;

					Ref = normal.unitVector().minus(p.direction.unitVector().divideScaler(2*cphi));

					//specular
					Ir += p.Ipr / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
					Ig += p.Ipg / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
					Ib += p.Ipb / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
				}
				else //case 3
				{
					//System.out.println("here2");
					Ref = p.direction.unitVector().flipDirection();

					//specular
					Ir += p.Ipr / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
					Ig += p.Ipg / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
					Ib += p.Ipb / dist * Ks * Math.pow((Ref.unitVector().dot(view)), n);
				}
			}

			//System.out.println("<" + Ir + ", " + Ig + ", " + Ib + ">");
			return new Vector3D(Ir, Ig, Ib);
		}*/

		public Vector3D isTranslucent()
		{
			boolean bColor;

			double Kdr = 0, Kdg = 0, Kdb = 0;

			if(intersection.X >= 0)
				bColor = true;
			else
				bColor = false;

			if(Math.abs(intersection.X%400) > 200)
				bColor = !bColor;
			if(Math.abs(intersection.Z%400) > 200)
				bColor = !bColor;

			if(bColor) //red
			{
				return null;
			}
			else //white
			{
				Kdr = transmitted * 1;
				Kdg = transmitted * 1;
				Kdb = transmitted * 1;

				return new Vector3D(Kdr, Kdg, Kdb);
			}
		}

		public double getTransmitted()
		{
			return transmitted;
		}
	}

	//and I said: "Let there be PointLight!"
    public class PointLight
    {
    	double distance;
    	float Ipr, Ipg, Ipb;
    	Vector3D direction;

    	public PointLight(double distance, double Ipr, double Ipg, double Ipb , Vector3D direction)
    	{
    		this.distance = distance;
    		this.Ipr = (float)Ipr;
    		this.Ipg = (float)Ipg;
    		this.Ipb = (float)Ipb;
    		this.direction = direction;
    	}

    	public PointLight(PointLight p, Vector3D direction)
    	{
    		this.distance = p.distance;
    		this.Ipr = p.Ipr;
    		this.Ipg = p.Ipg;
    		this.Ipb = p.Ipb;
    		this.direction = direction;
    	}

    	public void setDirection(double x, double y, double z)
    	{
    		direction.end.X = x;
    		direction.end.Y = y;
    		direction.end.Z = z;
    	}
    }
}