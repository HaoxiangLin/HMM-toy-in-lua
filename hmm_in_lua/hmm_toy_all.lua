--[[
  The MIT License
  
  Copyright (c) 2016, Chris Lin <chrislhx@gmail.com>
  
  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:
  
  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.
  
  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
]]--
 

-- configs 
local iter_times  = 500
local counts = {}
counts["ABBA"] = 10
counts["BAB"] = 20

local state = {"s", "t"}
local obs = {"A", "B"}

-- begin prob
local start_p = {}
start_p["s"] = 0.85
start_p["t"] = 0.15

local em_p = {}
em_p["s"] = {}
em_p["t"] = {}

em_p["s"]["A"] = 0.4
em_p["s"]["B"] = 0.6

em_p["t"]["A"] = 0.5
em_p["t"]["B"] = 0.5

local trans_p = {}
trans_p["s"] = {}
trans_p["t"] = {}

trans_p["s"]["t"] = 0.7
trans_p["t"]["s"] = 0.1 
trans_p["s"]["s"] = 0.3
trans_p["t"]["t"] = 0.9 

--------------funcitons 

function forward(target_seq,start_p,trans_p, em_p, alpha_mtx) 
-- alpha_mtx can be {}, when only use for calculate prob not use for bw
-- calcualte alpha in the BW algorithm
   local sum_pro = 0 
   local seq = {}
   for i = 1, #target_seq do 
      seq[#seq + 1] = target_seq:sub(i,i) 
   end  
   for k,_ in pairs(start_p) do 
       alpha_mtx[k] = {}
   end  
   for k,v in pairs(start_p) do  
      alpha_mtx[k][1] = start_p[k] * em_p[k][seq[1]] 
   end
   for i = 2, #target_seq do 
      for k,v in pairs(start_p) do  -- now
         local sum = 0
         for k1,v1 in pairs(start_p) do  -- before 
             sum = sum + em_p[k][seq[i]] * trans_p[k1][k] * alpha_mtx[k1][(i - 1)]
----            print (seq[i], k,k1,i,sum, em_p[k1][seq[i]] * trans_p[k][k1] * prob_mtx[k][(i - 1)], em_p[k1][seq[i]], trans_p[k][k1], prob_mtx[k][(i - 1)])
         end 
         alpha_mtx[k][i] = sum 
      end 
   end
   for k,v in pairs(alpha_mtx) do 
       for k1,v1 in ipairs(v) do 
--           print (k,k1,v1)
       end 
       sum_pro = sum_pro + alpha_mtx[k][#target_seq] 
   end  
   return sum_pro
end 
function viterbi(target_seq,start_p,trans_p,em_p)
   local max_likelihood = 0
   local max_likelihood_state_seq = {}
   local seq = {}
   for i = 1, #target_seq do 
      seq[#seq + 1] = target_seq:sub(i,i) 
   end  
   local prob_mtx = {}
   local back_track_mtx = {}
   for k,_ in pairs(start_p) do 
       prob_mtx[k] = {}
       back_track_mtx[k] = {}
   end  
   for k,v in pairs(start_p) do  
      prob_mtx[k][1] = start_p[k] * em_p[k][seq[1]] 
   end
 
   for i = 2, #target_seq do 
      for k,v in pairs(start_p) do  -- now
         local max_p_before = 0
         local max_p_state = ""
         for k1,_ in pairs(start_p) do  -- before 
             if trans_p[k1][k] * prob_mtx[k1][(i-1)] > max_p_before then 
                max_p_state = k1
                max_p_before =  trans_p[max_p_state][k] * prob_mtx[max_p_state][(i-1)]
             end
         end
         back_track_mtx[k][i] = max_p_state
         prob_mtx[k][i] = em_p[k][seq[i]] *  max_p_before 
  --           print (seq[i], k,k1,i,sum, em_p[k1][seq[i]] * trans_p[k][k1] * prob_mtx[k][(i - 1)], em_p[k1][seq[i]], trans_p[k][k1], prob_mtx[k][(i - 1)])
      end 
   end
   for k,v in pairs(prob_mtx) do 
     --  for k1,v1 in ipairs(v) do 
       if max_likelihood < prob_mtx[k][#target_seq] then 
          max_likelihood =  prob_mtx[k][#target_seq]   
          max_likelihood_state_seq[1] = k
       --    print (k,k1,math.log(v1)/math.log(2))
       end 
   end 
   for i = #target_seq, 1, -1  do 
--       print (max_likelihood_state_seq[#max_likelihood_state_seq], i,back_track_mtx[max_likelihood_state_seq[#max_likelihood_state_seq]][i]) 
       max_likelihood_state_seq[#max_likelihood_state_seq + 1] = back_track_mtx[max_likelihood_state_seq[#max_likelihood_state_seq]][i]
   end 
   return string.reverse(table.concat(max_likelihood_state_seq, "")),math.log(max_likelihood)/math.log(2)
end
function beta(target_seq,start_p,trans_p, em_p, len, beta_mtx)
--   print ("len", len)
   local seq = {}
   for i = 1, len do 
      seq[#seq + 1] = target_seq:sub(i,i) 
   end  
   for k,_ in pairs(start_p) do 
       beta_mtx[k] = {}
   end 
  -- beta should  be 1 at the end 
   for k,v in pairs(start_p) do  
      beta_mtx[k][len] = 1 
--      print (k, len, beta_mtx[k][len])
   end
   
   for i = len -1 , 1, -1 do 
      for k,v in pairs(start_p) do  -- now
         local sum = 0
         for k1,v1 in pairs(start_p) do  -- after 
            sum = sum + em_p[k1][seq[(i+1)]] * trans_p[k][k1] * beta_mtx[k1][(i+1)]
--            print (seq[(i+1)], k,k1,i,sum, em_p[k1][seq[i]] * trans_p[k][k1] * beta_mtx[k][(i - 1)], em_p[k1][seq[i]], trans_p[k][k1], prob_mtx[k][(i - 1)])
         end 
          beta_mtx[k][i] = sum 
--          print (k, i, beta_mtx[k][i])
      end 
   end
   for k,v in pairs(beta_mtx) do 
       for k1,v1 in ipairs(v) do 
--           print ("i",k,k1,v1)
       end 
   end  
end 
function gamma(target_seq,start_p,trans_p, em_p, alpha_mtx, beta_mtx, gamma_mtx, p_seq) 
-- calcualte alpha in the BW algorithm
   local sum_pro = 0 
   local seq = {}
   for i = 1, #target_seq do 
      seq[#seq + 1] = target_seq:sub(i,i) 
   end  
   for k,_ in pairs(start_p) do 
       gamma_mtx[k] = {}
       for k1,_ in pairs(start_p) do 
          gamma_mtx[k][k1] = {}
       end 
   end  
   for i = 1, #target_seq - 1 do -- can only calculate i+1 <= seq len 

      for k,v in pairs(start_p) do  -- now
         for k1,v1 in pairs(start_p) do  -- before 
             gamma_mtx[k1][k][i] = alpha_mtx[k1][i] * trans_p[k1][k] * em_p[k][seq[(i+1)]] * beta_mtx[k][(i+1)] / p_seq 

--            print  (k1, k ,i,gamma_mtx[k1][k][i], alpha_mtx[k1][i],trans_p[k1][k],em_p[k][seq[(i+1)]],beta_mtx[k][(i+1)], p_seq)
         end 
      end 
   end
--   for k,v in pairs(gamma_mtx) do 
--       for k1,v1 in ipairs(v) do 
--           for i = 1, #v1 do 
--              print (k1,k,v1[i])
--           end 
 --      end 
--       sum_pro = sum_pro + alpha_mtx[k][#target_seq] 
--   end  
   return sum_pro
end 
function delta(len, alpha_mtx, gamma_mtx, delta_mtx, p_seq)
   for k,_ in pairs(alpha_mtx) do 
       delta_mtx[k] = {}
       delta_mtx[k][len] = alpha_mtx[k][len]/p_seq
   end
   for i = 1,len-1 do 
      for k1,v1 in pairs(gamma_mtx) do  -- k1 before k after
         local sum = 0
         for k,v in pairs(v1) do 
--            print (k1,k,i,v[i])
            sum = sum + v[i]
         end 
         delta_mtx[k1][i]  = sum 
      end
   end  
--   for k,v in pairs(delta_mtx) do 
--      for k1,v1 in ipairs(v) do 
--         print (k, k1, v1)
--      end 
--   end 
end 

-- main ---------------------------------------

for iter = 1, iter_times do  
   local data = {}
   
   for k,v in pairs(counts) do 
      data[k] = {}
   --   print("process", k, string.len(k))
      data[k]["alpha"] = {} 
      data[k]["beta"] = {} 
      data[k]["gamma"] =  {}
      data[k]["delta"] =  {}
      local p_seq = forward(k,start_p,trans_p,em_p,data[k]["alpha"])
      beta(k,start_p,trans_p,em_p,string.len(k),data[k]["beta"])
      gamma(k, start_p,trans_p,em_p,data[k]["alpha"],data[k]["beta"],data[k]["gamma"], p_seq)
      delta(string.len(k), data[k]["alpha"], data[k]["gamma"], data[k]["delta"], p_seq)
   end
   
   --- cal start
   local new_start_counts = {}
   local new_start_sum = 0
   for _,v in pairs(state) do 
      for k1,v1 in pairs(counts) do
         new_start_counts[v] = (new_start_counts[v] or 0 ) + data[k1]["delta"][v][1] * v1 
      end 
      new_start_sum = new_start_sum + new_start_counts[v]
   end
   --print ("new start p:")
   for k,v in pairs(start_p) do
      start_p[k] = new_start_counts[k]/new_start_sum
   --   print(k,new_start_p[k], new_start_counts[k], new_start_sum)
   end
   
   --print ("new trans")
   for _,v1 in pairs(state) do
      local new_trans_counts = {}
      local new_trans_sum = 0
      for _,v2 in pairs(state) do
         
         for obs_seq, seq_count in pairs(counts) do
            local combine_p = 0
            for i = 1, #obs_seq  - 1 do 
               combine_p = combine_p + data[obs_seq]["gamma"][v1][v2][i] 
            end 
            new_trans_counts[v2] =  (new_trans_counts[v2] or 0) + combine_p * seq_count -- sum diff obs_seq
         end
         new_trans_sum =  new_trans_sum + new_trans_counts[v2] -- sum differ v1 -> v2
      end 
      for _,v2 in pairs(state) do
         trans_p[v1][v2] =  new_trans_counts[v2]/new_trans_sum
   --      print ( v1, "to", v2, new_trans_p[v1][v2], new_trans_counts[v2], new_trans_sum)
      end 
   end 
   
   --print ("new ems")
   for _,v1 in ipairs(state) do
      local new_em_counts = {}
      local new_em_sum = 0
      for _,symbol in ipairs(obs) do 
         for obs_seq, seq_count in pairs(counts) do
            local combine_p = 0
    --        print ("seq", obs_seq)
            for i = 1, #obs_seq do
               if obs_seq:sub(i,i) == symbol then
                  combine_p = combine_p + data[obs_seq]["delta"][v1][i]
      --           print (v1, symbol, i, data[obs_seq]["delta"][v1][i])
               end 
            end 
            new_em_counts[symbol] = ( new_em_counts[symbol] or 0) + combine_p * seq_count
         end
         new_em_sum = new_em_sum + new_em_counts[symbol]
      end 
      for _,symbol in ipairs(obs) do
         em_p[v1][symbol] = new_em_counts[symbol]/new_em_sum
   --      print (v1, symbol, new_em_p[v1][symbol], new_em_counts[symbol], new_em_sum)
      end 
   end
   
   local log_p_sum = 0
   for seq,seq_count in pairs(counts) do 
   --forward(k,start_p,trans_p,em_p,data[k]["alpha"])
      log_p_sum = log_p_sum + seq_count * math.log(forward(seq,start_p,trans_p,em_p,data[seq]["alpha"]))
   end 
   print ("iter", iter, "log likelihood",log_p_sum)
end 
 -- train finish 
--print (log_p_sum)
print(forward("ABBAAAABBB",start_p,trans_p,em_p,{}))
print(viterbi("ABBAAAABBB", start_p,trans_p,em_p))
--beta("ABBA",start_p,trans_p,em_p,4)
--beta("ABBA",start_p,trans_p,em_p,3)

