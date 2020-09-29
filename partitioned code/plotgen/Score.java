//Score Estimation Class
import java.util.*

public interface Tuple {
   int size();
}

public class Tuple2<T1,T2> implements Tuple {
   public final T1 item1;
   public final T2 item2;
 
   public Tuple2(final T1 item_1,final T2 item_2) {
      item1 = item_1;
      item2 = item_2;
   }

   @Override
   public int size() {
      return 2;
   }
}


public class Score(){
    public static HashMap<String,Float> getAllScores(int tp, int fp, int fn, int tn){
    	Map<String,Float> scoreMap = new HashMap<String,Float>;
        scoreMap["precision"] = (1.0 * tp) / (tp+fp)
    	scoreMap["recall"] = (1.0 * tp) / (tp+fn)
    	scoreMap["f1"] = (2.0 * scoreMap.get(precision) * scoreMap.get(recall)) / (scoreMap.get(precision) + scoreMap.get(recall)) 
        scoreMap["f01"] = (1.01 * scoreMap.get(precision) * scoreMap.get(recall)) / (0.01 * scoreMap.get(precision) + scoreMap.get(recall)) 
        return scoreMap
    }
    
    public static float getKendallTauScore(HashMap<Integer,HashMap<String,Set<Integer>>> truth, HashMap<Integer,HashMap<String,Set<Integer>>> predicted, ArrayList<String> scorestates){
        Set<Integer> allTimes = new HashSet<Integer>(truth.keySet());
        allTimes.addAll(predicted.keySet());
        float avgscore = 0.0; 
        for(String scorestate : scorestates) {
	   HashMap<Integer,Integer> node2TrueTime = new Map<Integer,Integer>();
           HashMap<Integer,Integer> node2MyTime = new Map<Integer,Integer>();  
           for(Integer time : truth.keySet()){
               for(Integer node : truth.get(time).get(scorestate)){
                   assert node2TrueTime.get(node).get(scorestate) != null;
		   node2TrueTime.put(node,time);
           }   
           for(Integer time : predicted.keySet()){
               for(Integer node : predicted.get(time).get(scorestate)){
                   assert node2MyTime.get(node).get(scorestate) != null;
		   node2MyTime.put(node,time);
           }
           int con, discon = 0, 0;  
           final int NOTSEEN = 1000000;
           Set<Integer> allNodes = new HashSet<Integer>(node2TrueTime.keySet());
           for(Integer node : allNodes){
	       if(node2MyTime.get(node) == null)
		  node2MyTime.put(node,NOTSEEN); 
               if(node2TrueTime.get(node) == null)
		  node2TrueTime.put(node,NOTSEEN) 
           }
           ArrayList<Integer> nodelist = new ArrayList<>(allNodes);           
           for(int i=0; i < nodelist.size(); i++){
	      int node1 = nodelist[i];
	      int true1, my1 = node2TrueTime.get(node1), node2MyTime.get(node1); 
	      for(int j=i+1; j < nodelist.size(); j++){
	          int node2 = nodelist[j];
                  int true2, my2 = node2TrueTime.get(node2), node2MyTime.get(node2); 
                  int signal = (true1 - true2) * (my1 - my2); 
                  if(signal >=0 ){
		     con += 1 ;
                  } else{
		     discon +=1 ;
		  }
              } 
           }    
           avgscore += 1.0 * (con - discon) / Math.sqrt((allNodes.size() - con) * (allNodes.size() - discon))  ;    
        }  
        return avgscore / scorestates.size();         
             
    }  
    
    public static float getMatchingScore(HashMap<Integer,HashMap<String,Set>> truth, HashMap<Integer,HashMap<String,Set>> predicted){
        return                  
    }
    
    public static HashMap<String,Float> estimateGraphScore(ArrayList<Tuple2<Integer,Integer>> truth, ArrayList<Tuple2<Integer,Integer>> predicted, String scorename="all"){
	assert ArrayList("all","f1","f01","precision","recall").exists(scorename);
       
        Set<Tuple2<Integer,Integer>> trueSet = new TreeSet<Tuple2<Integer,Integer>>();
        Set<Tuple2<Integer,Integer>> predictedSet = new TreeSet<Tuple2<Integer,Integer>>();
        Set<Integer> allNodes = new TreeSet<Integer>();
        
        Iterator itr = truth.iterator();
        while(itr.hasNext()) {
           Tuple2<Integer,Integer> trueItem = itr.next();
           trueSet.add(trueItem);
           allNodes.add(trueItem.item1)
           allNodes.add(trueItem.item2)
        }
        itr = predicted.iterator();
        while(itr.hasNext()) {
           Tuple2<Integer,Integer> preItem = itr.next();
           predictedSet.add(preItem);
           allNodes.add(preItem.item1)
           allNodes.add(preItem.item2)
        }
         
        int tp = trueSet.difference(predictedSet);
        int fp = trueSet.difference(predictedSet);
        int fn = trueSet.difference(predictedSet)l
	int tn = allNodes.length() * ( allNodes.length()-1 ) - tp - fp - fn;
	HashMap<String,Float> scores = this.getAllScores(int tp, int fp, int fn, int tn);
        return score
   }

}
