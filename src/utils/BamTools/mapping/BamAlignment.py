import sys
class SimplePropertyMapping:
    def __init__(self, mapstr):
        (self._dest_type, self._dest_name), (self._src_type, self._src_name) = map(lambda x: x.strip().split() , mapstr.split("==>")) 
    def name(self):
        return self._dest_name
    def class_def(self):
        return """
struct _{dest_name}_t {{
    _{dest_name}_t() : _ptr(NULL) {{}}
    void set(bam1_t* ptr) {{_ptr = ptr;}}
    void set(const _{dest_name}_t& that) {{_ptr = that._ptr;}}
    operator {dest_type}() const {{return ({dest_type})(_ptr == NULL?0:_ptr->{src_name});}}
    const {dest_type}& operator=(const {dest_type}& val) {{
        if(NULL != _ptr) _ptr->{src_name} = ({dest_type})val;
        return val;
    }}
    private:
        bam1_t* _ptr;
}} {dest_name};""".format(src_type = self._src_type,
                          dest_type = self._dest_type,
                          src_name = self._src_name,
                          dest_name = self._dest_name)
class ReadOnlyPropertyMapping:
    def __init__(self, mapstr):
        (self._dest_type, self._dest_name), (self._src_type, self._src_name) = map(lambda x: x.strip().split() , mapstr.split("==>")) 
    def name(self):
        return self._dest_name
    def class_def(self):
        return """
struct _{dest_name}_t {{
    _{dest_name}_t() : _ptr(NULL) {{}}
    void set(bam1_t* ptr) {{_ptr = ptr;}}
    void set(const _{dest_name}_t& that) {{_ptr = that._ptr;}}
    operator {dest_type}() const {{return ({dest_type})(_ptr == NULL?0:({src_name});}}
    private:
        bam1_t* _ptr;
}} {dest_name};""".format(src_type = self._src_type,
                          dest_type = self._dest_type,
                          src_name = self._src_name,
                          dest_name = self._dest_name)


props = {}
for line in file(sys.argv[1]):
    line = line.strip().split(":")
    if len(line) == 1: 
        ptype = "simple"
    else:
        ptype = line[1].strip()
    args = line[0]
    if ptype == "simple":
        prop = SimplePropertyMapping(args)
        props[prop.name()] = prop
    elif ptype == "readonly":
        prop = ReadOnlyPropertyMapping(args)
        props[prop.name()] = prop

for p in props:
    print props[p].class_def()

print "void setup(bam1_t* bam)"
print "{"
for p in props:
    print "    {name}.set(bam);".format(name = p)
print "}"

print "void setup(const BamAlignment& bam)"
print "{"
for p in props:
    print "    {name}.set(bam.{name});".format(name = p)
print "}"
