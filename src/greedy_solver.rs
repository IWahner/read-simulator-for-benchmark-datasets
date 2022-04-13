use rand::{thread_rng, Rng};
use std::collections::HashMap;

pub struct GreedySolver{
}

impl Default for GreedySolver {
    fn default() -> Self {
         Self::new()
      }
  }
 
impl GreedySolver{
 
     pub fn new() -> Self {
         GreedySolver {
         }
    }

    pub fn variantselection(&self, genloci: i32, numberofvariants: i32, wanted_freq: f32) -> HashMap<String, f32> {
        let mut rng = thread_rng();
        let mut sol: HashMap<String,f32> = HashMap::new();
        //too avoid picking too much variants
        let mut break_condition = genloci as f32 * wanted_freq;
        //modify the freq for the number of variants
        let mod_freq = wanted_freq * (genloci as f32/numberofvariants as f32);
        // to get frequencies till 0,001 %
        let cut_off = mod_freq * 100000.0;
        println!("cutoff: {}", cut_off);
        for i in 0..numberofvariants{
            if break_condition == 0.0 {
                break;
            }
            let num: i32 = rng.gen_range(0..=100000);
            //println!("num: {}", num);
            if num <= cut_off as i32 {
                sol.insert("x".to_owned()+&i.to_string(), 1.0);
                break_condition = break_condition - 1.0;
            }
            else {
                sol.insert("x".to_owned()+&i.to_string(), 0.0);
            }
        }
        let mut sum_of_ones = 0.0;
        for key in sol.values(){
            sum_of_ones = sum_of_ones + key;
        }
        println!("Sum_of_ones: {}", sum_of_ones);
        println!("Länge: {}", sol.len());
        println!("Anzahl an errechnete Summe: {}", sum_of_ones/genloci as f32);
        sol
    }
}