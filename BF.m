Kv = [0.001:0.001:0.01, 0.02:0.01:1, 1.02:0.02:10];

Ti = Kv; Kp = Kv;

tmin = 1000;

kb= 0;
tb =0;
j = 0;
for i = 1: length(Kv)
for k = 1: length(Kv)
    Cc = -Kp(i)*(Ti(k)*s+1)/Ti(k)*s;
    
    o = stepinfo(feedback(Cc*Gp,1)).Overshoot;
    t = stepinfo(feedback(Cc*Gp,1)).SettlingTime;
    
    if o <= 1 && t <= tmin
        kb= Kp(i);
        tb =Ti(k);
    
    end
    
end
j = j+1;
disp(j)
end
