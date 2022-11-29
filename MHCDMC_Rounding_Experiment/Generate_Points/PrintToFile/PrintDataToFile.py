import os

def print_data_to_file(User, Server, fname):
    # 参数：
    # User      list    [[x1, x2, BR1], ..., [xn, yn, BRn]]
    # Server    list    [[x1, x2, BW1], ..., [xn, yn, BWm]]
    # fname         str     文件的绝对路径
    n = len(User)
    m = len(Server)
    with open(fname, 'a') as f:
        s = "%d %d\n" % (m, n)
        f.write(s)
        for a in Server:
            s = "%9f %9f %d\n" % (a[0], a[1], a[2])
            f.write(s)
        for u in User:
            s = "%f %f %d\n" % (u[0], u[1], u[2])
            f.write(s)
    pass