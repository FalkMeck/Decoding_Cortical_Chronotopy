noise_50 = c(73,51,26,75,98,41,44,67,20,3,65,31,52,22,42,66,25,61,50,87)
isi_dur = c(3250,3450,3650,3850,4050,4250,4450,4650,4850,5050,5450,5650,5850,6050,6250,6450,6650,6850,7050,7250)

#find_5 = FALSE
#while (!find_5) {
#  noise_5 = sample(noise_50,5)
#  isi_5 = sample(isi_dur,5)
#  pause_mean = mean(noise_5+isi_5)
#  find_5 = pause_mean == (7400-2100)
#}

noise_5 = c(44,61,3,75,67) # maybe have to change one of the to a millisecond slower number
isi_5 = c(3850,4050,5450,7050,5850)

noise_full = c()
isi_full = c()
for (k in 1:4) {
while_dur = 0
while (while_dur != 1) {
  noise_5_shuffle = sample(noise_5)
  isi_5_shuffle = sample(isi_5)
  cumsum = 0
  unequal_counter = 0
  for (i in 1:length(isi_5_shuffle)) {
    pause_total = noise_5_shuffle[i] + isi_5_shuffle[i]+2100
    if (cumsum/1000 != floor(cumsum/1000)){
      unequal_counter = unequal_counter +1
    }
    cumsum = cumsum + 162*700
  }
  if (unequal_counter == length(isi_5_shuffle)) {
    while_dur = 1
  }
}
noise_full = cbind(noise_full, noise_5_shuffle)
isi_full = cbind(isi_full, isi_5_shuffle)
}
