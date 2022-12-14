function v_ol = overlap(s1,s2,nw,trSpan)
    v_ol = exp(-min([abs(s1-s2),trSpan-abs(s1-s2)])^2/(4*nw^2));
end