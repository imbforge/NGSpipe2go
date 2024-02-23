class Ngspipe2goWrongTypeException extends Exception {
  Ngspipe2goWrongTypeException(String message) {
    super(message)
  }
}

Boolean validate_schema(Class Params, Map params, String module) {
  try {
    // validate parameter types against the schema
    p = Params.newInstance(params)
    params.each{ k, v -> 
      if(p[k].getClass() != params[k].getClass()) {
        String message = "param ${k} is ${params[k].getClass()} instead of ${p[k].getClass()}"
        throw new Ngspipe2goWrongTypeException(message)
      }
    }
    // validate presence of mandatory parameters
    assert true == !!p
  } catch(Ngspipe2goWrongTypeException e) {
    throw new RuntimeException("invalid parameter types in module ${module}\n${e}")
  } catch(AssertionError e) {
    throw new RuntimeException("mandatory arguments missing in module ${module}")
  }
  return true
}
