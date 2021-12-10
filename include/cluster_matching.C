//function to order clust_x_i and clust_y_i indices by largest to smallest in ADC value
void cluster_matching(vector<float> &clust_x, vector<float> &clust_y, vector<int> &clust_x_i, vector<int> &clust_y_i){


  for(int i = 0; i < clust_x.size(); i++){
    double largest = clust_x[i];
    int largest_i = i;

    for(int j = i+1; j < clust_x.size(); j++){
      if(clust_x[j] > largest){
  largest = clust_x[j];
  largest_i = j;
      }
    }

    clust_x[largest_i] = clust_x[i];
    clust_x[i] = largest;
    
    int temp_i = clust_x_i[largest_i];
    clust_x_i[largest_i] = clust_x_i[i];
    clust_x_i[i] = temp_i;
  }

  for(int i = 0; i < clust_y.size(); i++){
    double largest = clust_y[i];
    int largest_i = i;
    
    for(int j = i+1; j < clust_y.size(); j++){
      if(clust_y[j] > largest){
        largest = clust_y[j];
        largest_i = j;
      }
    }
    
    clust_y[largest_i] = clust_y[i];
    clust_y[i] = largest;

    int temp_i = clust_y_i[largest_i];
    clust_y_i[largest_i] = clust_y_i[i];
    clust_y_i[i] = temp_i;
  }


}
