function d = point_to_line(a, b, c)

ab=sqrt((a(1,1)-b(1,1))^2+(a(1,2)-b(1,2))^2);
ac=sqrt((a(1,1)-c(1,1))^2+(a(1,2)-c(1,2))^2);
bc=sqrt((c(1,1)-b(1,1))^2+(c(1,2)-b(1,2))^2);
cross=(b(1)-c(1))*(a(1)-c(1))+(b(2)-c(2))*(a(2)-c(2));
cos_theta=(ab^2+bc^2-ac^2)/(2*ab*bc);
if cross<0 || cross>bc
    d=min(sqrt((a(1)-b(1))^2+(a(2)-b(2))^2),sqrt((a(1)-c(1))^2+(a(2)-c(2))^2));
else 
    d=ab*sqrt(1-cos_theta*cos_theta);
end

