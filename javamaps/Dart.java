import java.util.ArrayList;
import java.util.LinkedList;

class Dart{
	
    int dim;
    Dart[] pointers;
    Vector[] attributes;

    public static int modulo(int a, int b){
        return (((a % b) + b) % b);
    }

    public static void main(String[] args){
        Dart[] darts = new Dart[12];
        for(int i=0; i<darts.length; i++){
            darts[i] = new Dart(2);
            darts[i].setDAttribute(0, new Vector(i));
        }

        // Test 2: One tetrahedron
        makeTriangle(darts[0], darts[1], darts[2]);
        makeTriangle(darts[3], darts[4], darts[5]);
        makeTriangle(darts[6], darts[7], darts[8]);
        makeTriangle(darts[9], darts[10], darts[11]);

        dSew(2, darts[5], darts[6]);
        dSew(2, darts[0], darts[4]);
        dSew(2, darts[2], darts[7]);
        dSew(2, darts[3], darts[9]);
        dSew(2, darts[8], darts[10]);
        dSew(2, darts[1], darts[11]);

        ArrayList<Dart> cell0 = darts[5].getDCell(0);
        for(int i=0; i<cell0.size();i++){
            System.out.println(cell0.get(i).getDAttribute(0).dim);
        }

        System.out.println();

        cell0 = darts[5].getDCell(1);
        for(int i=0; i<cell0.size();i++){
            System.out.println(cell0.get(i).getDAttribute(0).dim);
        }

        System.out.println();

        cell0 = darts[5].getDCell(2);
        for(int i=0; i<cell0.size();i++){
            System.out.println(cell0.get(i).getDAttribute(0).dim);
        }

        System.out.println();

        cell0 = darts[5].getDCell(3);
        for(int i=0; i<cell0.size();i++){
            System.out.println(cell0.get(i).getDAttribute(0).dim);
        }
        System.out.println();
        for(int i=0; i<darts.length; i++){
            System.out.println(darts[i]);
        }

        // Test 1: Planar graph
        // darts[0].setDPointer(0, darts[2]);
        // darts[0].setDPointer(1, darts[1]);
        // darts[1].setDPointer(0, darts[0]);
        // darts[1].setDPointer(1, darts[2]);
        // darts[2].setDPointer(0, darts[1]);
        // darts[2].setDPointer(1, darts[0]);
        // darts[2].setDPointer(2, darts[3]);
        // darts[3].setDPointer(0, darts[4]);
        // darts[3].setDPointer(1, darts[6]);
        // darts[3].setDPointer(2, darts[2]);
        // darts[4].setDPointer(0, darts[5]);
        // darts[4].setDPointer(1, darts[3]);
        // darts[5].setDPointer(0, darts[6]);
        // darts[5].setDPointer(1, darts[4]);
        // darts[6].setDPointer(0, darts[3]);
        // darts[6].setDPointer(1, darts[5]);

        // ArrayList<Dart> cell1 = darts[1].getDCell(1);
        // for(int i=0; i<cell1.size(); i++){
        //     System.out.println(cell1.get(i).getDAttribute(0).dim);
        // }

        // ArrayList<Dart> cell2 = darts[1].getDCell(2);
        // for(int i=0; i<cell2.size(); i++){
        //     System.out.println(cell2.get(i).getDAttribute(0).dim);
        // }

        // ArrayList<Dart> cell3 = darts[1].getDCell(3);
        // for(int i=0; i<cell3.size(); i++){
        //     System.out.println(cell3.get(i).getDAttribute(0).dim);
        // }
        
        // ArrayList<Dart> altcell1 = darts[2].getDCell(1);
        // for(int i=0; i<altcell1.size(); i++){
        //     System.out.println(altcell1.get(i).getDAttribute(0).dim);
        // }

        // ArrayList<Dart> altcell2 = darts[4].getDCell(2);
        // for(int i=0; i<altcell2.size(); i++){
        //     System.out.println(altcell2.get(i).getDAttribute(0).dim);
        // }
    }

    public Dart(int d){
        dim = d;
        pointers = new Dart[d+1];
        attributes = new Vector[d+1];
    }

    // Attaches the d-face containing dart1 and the d-face containing dart2
    // by connecting two of their d-1 faces (making them d-point to each other)
    public static boolean dSew(int d, Dart dart1, Dart dart2){
        ArrayList<Dart> cell1 = dart1.getDCell(d-1);
        ArrayList<Dart> cell2 = dart2.getDCell(d-1);
        if(cell1.size() != cell2.size()){
            return false;
        }
        int size = cell1.size();
        int cell2_index = 0;
        for(int i=0; i<size; i++){
            cell1.get(i).setDPointer(d, cell2.get(cell2_index));
            cell2.get(cell2_index).setDPointer(d, cell1.get(i));
            cell2_index -= 1;
            if (cell2_index < 0){
                cell2_index = modulo(cell2_index, size);
            }
        }
        return true;
    }

    public static void makeTriangle(Dart dart1, Dart dart2, Dart dart3){
        dart1.setDPointer(0, dart3);
        dart1.setDPointer(1, dart2);
        dart2.setDPointer(0, dart1);
        dart2.setDPointer(1, dart3);
        dart3.setDPointer(0, dart2);
        dart3.setDPointer(1, dart1);
    }

    // Sews a set of at least three darts into a ccw polygon
    public static void makePolygon(Dart[] darts){
        for(int i=1; i<darts.length-1; i++){
            darts[i].setDPointer(0, darts[i-1]);
            darts[i].setDPointer(1, darts[i+1]);
        }
        darts[0].setDPointer(0, darts[darts.length-1]);
        darts[0].setDPointer(1, darts[1]);
        darts[darts.length-1].setDPointer(0, darts[darts.length-2]);
        darts[darts.length-1].setDPointer(1, darts[0]);
    }

    public void setDPointer(int d, Dart dart1){
        pointers[d] = dart1;
    }

    public Dart getDPointer(int d){
        return pointers[d];
    }

    public void setDAttribute(int d, Vector v){
        attributes[d] = v;
    }

    public Vector getDAttribute(int d){
        return attributes[d];
    }

    public String toString(){
        String rtn = "";
        for(int i=0; i<=dim; i++){
            rtn = rtn + "b" + i + ": " + (pointers[i]==null?"null":pointers[i].getDAttribute(0).dim) + "\n";
        }
        return rtn;
    }

    // Returns all the darts belonging to the same d-cell as this one
    public ArrayList<Dart> getDCell(int d){
        ArrayList<Dart> visited = new ArrayList<Dart>();
        LinkedList<Dart> to_visit = new LinkedList<Dart>();
        to_visit.add(this);
        visited.add(this);
        
        if(d==0){
            // Use breadth-first search to find the orbit of this dart under
            // 1) composited Bi Bj where 1<=i<j<=d
            // 2) the inverse of these functions - Bj Bi when neither i nor j is 1 and replacing the 1 with 0 otherwise
            while (to_visit.size() != 0){
                Dart curr = to_visit.poll();
                for(int i=1; i<=dim; i++){
                    for(int j=i+1; j<=dim; j++){
                        Dart newdart1 = curr.getDPointer(j);
                        if (newdart1 != null){
                            newdart1 = newdart1.getDPointer(i);
                            if (newdart1 != null && !visited.contains(newdart1)){
                                to_visit.add(newdart1);
                                visited.add(newdart1);
                            }
                        }
                        
                        Dart newdart2; 
                        if (i==1){
                            newdart2 = curr.getDPointer(0);
                        } else {
                            newdart2 = curr.getDPointer(i);
                        }
                        if (newdart2 != null){
                            newdart2 = newdart2.getDPointer(j);
                            if(newdart2 != null && !visited.contains(newdart2)){
                                to_visit.add(newdart2);
                                visited.add(newdart2);
                            }
                        }
                    }
                }
            }
        } else {
            // Use breadth-first search to find the orbit of this dart under
            // Bi for 1<=i<=dim; i!=d
            while(to_visit.size() != 0){
                Dart curr = to_visit.poll();
                for(int i=1; i<=dim; i++){
                    if (i!=d){
                        Dart newdart = curr.getDPointer(i);
                        if (newdart != null && !visited.contains(newdart)){
                            to_visit.add(newdart);
                            visited.add(newdart);
                        }
                    }
                }
            }
        }
        return visited;
    }

}