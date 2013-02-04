l = []
l.append((0.03,0.03))
l.append((0.03,0.97))
l.append((0.50,0.97))
l.append((0.50,0.50))
l.append((0.53,0.50))
l.append((0.53,0.97))
l.append((0.97,0.97))
l.append((0.97,0.03))
ow = (230,230)
#if len(sys.argv) > 1:
#    ow[0] = int(sys.argv[0])
#ow[1] = int(sys.argv[1])

mask = ""
mask += str(ow[0]) + "\n" + str(ow[1]) + "\n230.0\n230.0\n10.0\n0.1\n450.0\n350.0\n300.0\n2.0\n0.1064\n0.0000001072\n0.00001172\n350.0\n"
for x in range(ow[0]+1):
    for y in range(ow[1]):
        # determine if a point is inside a given polygon or not
        # Polygon is a list of (x,y) pairs.
        n = len(l)
        
        inside =False
        p1x, p1y = l[0]
        p1x *= ow[0];
        p1y *= ow[1];
        for i in range(n+1):
            p2x,p2y = l[i % n]
            p2x *= ow[0];
            p2y *= ow[1];
            if y > min(p1y,p2y):
                if y <= max(p1y,p2y):
                    if x <= max(p1x,p2x):
                        if p1y != p2y:
                            xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x,p1y = p2x,p2y

        if inside == True:
            mask+="1"
        else:
            mask+="0"
    mask+="\n"
print(mask)

        

        
