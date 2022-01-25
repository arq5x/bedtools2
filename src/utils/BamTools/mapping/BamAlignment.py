import sys
class SimplePropertyMapping:
    def __init__(self, mapstr):
        (self._dest_type, self._dest_name), (self._src_type, self._src_name) = map(lambda x: x.strip().split() , mapstr.split("==>")) 
    def name(self):
        return self._dest_name
    def class_def(self):
        return """
struct _{dest_name}_t {{
    operator {dest_type}() const {{return ({dest_type})(_ptr()->{src_name});}}
    const {dest_type}& operator=(const {dest_type}& val) {{
        _ptr()->{src_name} = ({dest_type})val;
        return val;
    }}
    private:
        bam1_t* _ptr() const {{
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->{dest_name})))->HtsObj2();
        }}
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
    operator {dest_type}() const {{return ({dest_type})({src_name});}}
    private:
        bam1_t* _ptr() const {{
           return ((BamAlignment*)(((uintptr_t)this) - ((uintptr_t)&((BamAlignment*)NULL)->{dest_name})))->HtsObj2();
        }}
}} {dest_name};""".format(src_type = self._src_type,
                          dest_type = self._dest_type,
                          src_name = self._src_name,
                          dest_name = self._dest_name)


props = {}
for line in open(sys.argv[1]):
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
    print(props[p].class_def())
