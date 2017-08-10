import java.io.*;
import java.lang.*;
import java.util.*;
import java.awt.*;
import java.math.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.FontMetrics;

public class annotate_tree_with_dup_modes{
String tree;
int n_node,n_leaf,xdim,ydim,font_size;
Hashtable <Integer,String> leaf1;
Hashtable <String,Integer> leaf2;
Hashtable <Integer,node> node_list;
Vector <branch> branch_list;
Hashtable <String, String>dup_pair;
Vector <String>family_pair;
double div_line=0.75;
public void read_synteny(String fpath)
{
try{
  FileInputStream fstream = new FileInputStream(fpath);
  DataInputStream in = new DataInputStream(fstream);
  BufferedReader br = new BufferedReader(new InputStreamReader(in));
  dup_pair=new Hashtable<String, String>();
  String strLine;
  strLine = br.readLine();
  while ((strLine = br.readLine()) != null)
  {
  if(strLine.length()==0)
  {
  continue;
  }
  if(strLine.indexOf("#epochs")>-1)
  {
  //System.out.println(strLine);
  break;
  }
  String[] ss=strLine.split("\t");
  String temp=ss[0]+","+ss[2];
  dup_pair.put(temp,ss[5]);
  }
  }
catch (Exception e)
   {
      System.err.println("Reading collinearity file error: " + e.getMessage());
   }
}
public void read_tree(String fpath)
{
try{
       FileInputStream fstream = new FileInputStream(fpath);
       DataInputStream in = new DataInputStream(fstream);
       BufferedReader br = new BufferedReader(new InputStreamReader(in));
       String strLine;
       String temp="";
       while ((strLine = br.readLine()) != null)
       {
         strLine.trim();
         if(strLine.length()>0)
         {
           temp=temp+strLine;
         }
       }
       tree=temp.substring(1,temp.length()-2);
////////////////////////////////////////////
  }
  catch (Exception e)
    {  
      System.err.println("Reading tree file error: " + e.getMessage());
    }
}
public rt recursion(int start, int end,double xpos)
{
   int[] brat={-1,-1,-1,-1,-1,-1};
   int i,j,num,tag;
   j=num=tag=0;
   for(i=start;i<=end;i++)
   {
     if(tree.charAt(i)=='(')
      {
       if(tag==1)
        {
        j++;
        }
        if(j==0)
        {
        brat[2*num]=i;
        }
        tag=1;
      }
      if(tree.charAt(i)==')')
      {
       if(tag==-1)
        {
        j--;
        }
        if(j==0)
        {
        brat[2*num+1]=i;
        num++;
        }
        tag=-1;
      }
   }
   String temp_tree="";
   if(num>0)
   {
     if(brat[0]>start)
       {
         temp_tree+=tree.substring(start,brat[0]);
       }
     for(i=0;i<num;i++)
       {
         if(i>0)
           {
           temp_tree+=tree.substring(brat[2*i-1]+1,brat[2*i]);
           }
         temp_tree+="xxxxxxxx";
       }
     temp_tree+=tree.substring(brat[2*num-1]+1,end+1);
   }
   else
   {
     temp_tree=tree.substring(start,end+1);
   }
   String[] ss=temp_tree.split(",");
   int[] next_node=new int[ss.length];
   double ypos=0.0;
   i=j=0;
   for(i=0;i<ss.length;i++)
   {
   String[] tt=ss[i].split(":");
   if(tt[0].indexOf("xxxxxxxx")>-1)
     {
       rt temp_rt=recursion(brat[2*j]+1,brat[2*j+1]-1,xpos+Double.parseDouble(tt[1]));
       next_node[i]=temp_rt.node;
       ypos+=temp_rt.ypos;
       j++;
     }
   else
     {
       int index=leaf2.get(tt[0]);
       ypos+=node_list.get(index).ypos;
       node_list.get(index).xpos=xpos+Double.parseDouble(tt[1]);
       next_node[i]=index;
     }
   }
   for(i=0;i<ss.length;i++)
   {
     branch temp_branch=new branch(n_node,next_node[i]);
     branch_list.add(temp_branch);
   }
   double temp_y=ypos/(double)ss.length;
   node temp_node=new node(xpos,temp_y);
   node_list.put(n_node,temp_node);
   rt value=new rt(n_node,temp_y);
   n_node++;
   return value;
}
public void processtree()
{
  char prev=' ';
  String id="";
  int i,j,start1;
  i=j=start1=0;
  leaf1=new Hashtable<Integer,String>();
  leaf2=new Hashtable<String,Integer>();
  node_list=new Hashtable<Integer,node>();
  branch_list=new Vector<branch>();
  node temp_node;
  while(i<tree.length())
  {
  if(tree.charAt(i)!='('&&tree.charAt(i)!=')')
  {
     if(prev=='('||prev==',')
     {
     id=id+tree.charAt(i);
     start1=1;
     }
     else
     {
       if(start1==1)
       {
         if(tree.charAt(i)==':')
         {
         leaf1.put(j,id);
         leaf2.put(id,j);
         temp_node=new node(-1.0,(double)j);
         node_list.put(j,temp_node);
         id="";
         start1=0;
         j++;
         }
         else
         {
         id=id+tree.charAt(i);
         }
       }
     }
  }
  prev=tree.charAt(i);
  i++;
  }
  n_node=n_leaf=j;
  recursion(0,tree.length()-1,0.0);
}
public void comp_pairs()
{
  family_pair=new Vector<String>();
  int i,k;
  String temp;
  for(i=0;i<leaf1.size()-1;i++)
  {
  for(k=i+1;k<leaf1.size();k++)
  {
  temp=leaf1.get(i)+","+leaf1.get(k);
  if(dup_pair.containsKey(temp))
  {
  family_pair.add(temp);
  }
  temp=leaf1.get(k)+","+leaf1.get(i);
  if(dup_pair.containsKey(temp))
  {
  family_pair.add(temp);
  }
  }
  }
}
public void paint (Graphics g) {
    Font font1 = new Font("Helvetica", Font.PLAIN,  font_size);
    g.setFont(font1);
    FontMetrics fm=g.getFontMetrics();
    g.setColor(Color.white);
    g.fillRect(0,0,xdim,ydim);
    g.setColor(Color.black);
    Hashtable <Integer,Integer> leaf_wid;
    leaf_wid=new Hashtable <Integer, Integer>();
    int i,j;
    double max_len=0;
    int max_wid=0;
    for(i=0;i<n_leaf;i++)
    {
    if(node_list.get(i).xpos>max_len)
    max_len=node_list.get(i).xpos;
    int temp=fm.stringWidth(leaf1.get(i));
    leaf_wid.put(i,temp);
    if(temp>max_wid)
    {
    max_wid=temp;
    }
    }
    int xmargin=(int)((float)xdim*0.05);
    int ymargin=(int)((float)ydim*0.05);
    int wordwidth;
    int charHeight=fm.getHeight();
    int xaxis=(int)((float)xdim*div_line-xmargin-(float)(max_wid)*1.05);
    //System.out.println(xaxis);
    int yaxis=(int)((float)ydim*0.9);
    g.drawString("Segmental",(int)((float)xdim*0.11),(int)((float)ydim*0.025+(float)charHeight/2));
    wordwidth=fm.stringWidth("Segmental");
    g.drawString("Tandem",(int)((float)xdim*0.11+wordwidth+(float)xdim*0.07),(int)((float)ydim*0.025+(float)charHeight/2));
    wordwidth=fm.stringWidth("Segmental tandem");
    g.drawString("Proximal",(int)((float)xdim*0.11+wordwidth+(float)xdim*0.14),(int)((float)ydim*0.025+(float)charHeight/2));
    wordwidth=fm.stringWidth("Segmental tandem proximal");
    g.drawString("Transposed",(int)((float)xdim*0.11+wordwidth+(float)xdim*0.21),(int)((float)ydim*0.025+(float)charHeight/2));
    //wordwidth=fm.stringWidth("Segmental tandem proximal transposed");
    //g.drawString("Ancestral locus",(int)((float)xdim*0.11+wordwidth+(float)xdim*0.28),(int)((float)ydim*0.025+(float)charHeight/2));
    g.setColor(Color.red);
    g.drawLine((int)((float)xdim*0.05),(int)((float)ydim*0.025),(int)((float)xdim*0.1),(int)((float)ydim*0.025));
    g.setColor(Color.orange);
    wordwidth=fm.stringWidth("Segmental");
    g.drawLine((int)((float)xdim*0.05+wordwidth+(float)xdim*0.07),(int)((float)ydim*0.025),(int)((float)xdim*0.1+wordwidth+(float)xdim*0.07),(int)((float)ydim*0.025));
    g.setColor(Color.green);
    wordwidth=fm.stringWidth("Segmental tandem");
    g.drawLine((int)((float)xdim*0.05+wordwidth+(float)xdim*0.14),(int)((float)ydim*0.025),(int)((float)xdim*0.1+wordwidth+(float)xdim*0.14),(int)((float)ydim*0.025));
g.setColor(Color.blue);
    wordwidth=fm.stringWidth("Segmental tandem proximal");
    g.drawLine((int)((float)xdim*0.05+wordwidth+(float)xdim*0.21),(int)((float)ydim*0.025),(int)((float)xdim*0.1+wordwidth+(float)xdim*0.21),(int)((float)ydim*0.025));
    g.drawLine((int)((float)xdim*0.1+wordwidth+(float)xdim*0.21),(int)((float)ydim*0.025),(int)((float)xdim*0.1+wordwidth+(float)xdim*0.21-(float)charHeight/2),(int)((float)ydim*0.025-(float)charHeight/2));
    g.drawLine((int)((float)xdim*0.1+wordwidth+(float)xdim*0.21),(int)((float)ydim*0.025),(int)((float)xdim*0.1+wordwidth+(float)xdim*0.21-(float)charHeight/2),(int)((float)ydim*0.025+(float)charHeight/2));
    //wordwidth=fm.stringWidth("Segmental tandem proximal transposed");
    //g.drawLine((int)((float)xdim*0.08+wordwidth+(float)xdim*0.28),(int)((float)ydim*0.025-(float)charHeight/2),(int)((float)xdim*0.08+wordwidth+(float)xdim*0.28),(int)((float)ydim*0.025+(float)charHeight/2));
    g.setColor(Color.black);
    for(i=0;i<n_leaf;i++)
    {
    g.drawString(leaf1.get(i),(int)((double)xaxis*node_list.get(i).xpos/max_len+xmargin),(int)((double)yaxis*(node_list.get(i).ypos+1.0)/(double)(n_leaf+1)+(double)ymargin+(double)charHeight/2.0));
    }    int x1,y1,x2,y2;
    for(i=0;i<branch_list.size();i++)
    {    
    x1=(int)((double)xaxis*node_list.get(branch_list.get(i).node1).xpos/max_len+xmargin);
    y1=(int)((double)yaxis*(node_list.get(branch_list.get(i).node1).ypos+1.0)/(double)(n_leaf+1))+ymargin;    x2=(int)((double)xaxis*node_list.get(branch_list.get(i).node2).xpos/max_len+xmargin);
    y2=(int)((double)yaxis*(node_list.get(branch_list.get(i).node2).ypos+1.0)/(double)(n_leaf+1))+ymargin;    g.drawLine(x1,y1,x1,y2);
    g.drawLine(x1,y2,x2,y2);    
    }
////////////////////////////////////////////////////////////
    //g.setColor(Color.red);
    double ny1,ny2,nx1,nx2;
    String []ss0;
    int id1,id2;
    for(i=0;i<family_pair.size();i++)
    {
    ss0=family_pair.get(i).split(",");
    id1=leaf2.get(ss0[0]);
    id2=leaf2.get(ss0[1]);
    nx1=(double)xaxis*node_list.get(id1).xpos/max_len+xmargin+leaf_wid.get(id1)+2;
    nx2=(double)xaxis*node_list.get(id2).xpos/max_len+xmargin+leaf_wid.get(id2)+2;
    ny1=(double)yaxis*(node_list.get(id1).ypos+1.0)/(double)(n_leaf+1)+(double)ymargin;
    ny2=(double)yaxis*(node_list.get(id2).ypos+1.0)/(double)(n_leaf+1)+(double)ymargin;
    String col=dup_pair.get(family_pair.get(i));
    if(col.compareTo("segmental")==0)
    g.setColor(Color.red);
     if(col.compareTo("tandem")==0)
    g.setColor(Color.orange);
if(col.compareTo("proximal")==0)
    g.setColor(Color.green);
if(col.compareTo("transposed")==0)
{
    g.setColor(Color.blue); 
//g.drawLine((int)(nx2),(int)(ny2-(float)charHeight/2),(int)(nx2),(int)(ny2+(float)charHeight/2));
g.drawLine((int)(nx1),(int)(ny1),(int)(nx1+(float)charHeight/2),(int)(ny1-(float)charHeight/2));
g.drawLine((int)(nx1),(int)(ny1),(int)(nx1+(float)charHeight/2),(int)(ny1+(float)charHeight/2));
} 
  drawBeizer(g,nx1, ny1,nx2,ny2);
    }
////////////////////////////////////////////////////////////   
}
public void drawBeizer(Graphics g,double x1,double y1,double x2,double y2)
   {
   double[] GX = new double[4];
   double[] GY = new double[4];
   GX[0]=x1;
   GX[1]=x1+(double)(xdim*(1-div_line));
   GX[2]=x2+(double)(xdim*(1-div_line));
   GX[3]=x2;
   GY[0]=y1;
   GY[1]=y1; 
   GY[2]=y2;
   GY[3]=y2;
   Cubic xSpline = new Cubic(Cubic.BEZIER, GX);
   Cubic ySpline = new Cubic(Cubic.BEZIER, GY);
   for (double t = 0 ; t < 1 ; t += 0.01)
   g.drawLine((int)xSpline.eval(t),(int)ySpline.eval(t),(int)xSpline.eval(t+0.01),(int)ySpline.eval(t+0.01));
   }
public static void main(String args[])
{
if(args.length<6)
{
System.out.println("Usage: java annotate_tree_with_dup_modes -t tree_file -s duplication_file(the output of detect_dup_modes_for_a_family.pl) -o output_PNG_file");
System.out.println("optional:-x plot_width -y plot height -f font_size");
System.exit(1);
}
HashMap<String,String> option = new HashMap<String, String>();
int i;
for(i=0;i<args.length/2;i++)
{
option.put(args[2*i],args[2*i+1]);
}
if(!option.containsKey("-t")||!option.containsKey("-s")||!option.containsKey("-o"))
{
System.out.println("Usage: java annotate_tree_with_dup_modes -t tree_file -s duplication_file(the output of detect_dup_modes_for_a_family.pl) -o output_PNG_file");
System.out.println("optional:-x plot_width -y plot height -f font_size");
System.exit(1);
}

annotate_tree_with_dup_modes proc=new annotate_tree_with_dup_modes();
proc.xdim=800;
proc.ydim=800;
proc.font_size=12;
if(option.containsKey("-x"))
proc.xdim=Integer.parseInt(option.get("-x"));
if(option.containsKey("-y"))
proc.ydim=Integer.parseInt(option.get("-y"));
if(option.containsKey("-f"))
proc.font_size=Integer.parseInt(option.get("-f"));
proc.read_synteny(option.get("-s"));
proc.read_tree(option.get("-t"));
proc.processtree();
proc.comp_pairs();
try  {
      BufferedImage bi = new BufferedImage(proc.xdim, proc.ydim, BufferedImage.TYPE_INT_ARGB);
      Graphics2D ig2 = bi.createGraphics();
      proc.paint(ig2);
      ImageIO.write(bi, "PNG", new File(option.get("-o")));
      } 
      catch (IOException ie)
      {ie.printStackTrace();}
}
}

