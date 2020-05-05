var resourcer_resource = function(type, name, params, credentials) {

  //
  // Resource factory functions to be reused
  //
  var toGridfsResource = function(name, params, credentials) {
      return {
          name: name,
          url: "gridfs://" + params.host + ":" + params.port + "/" + params.db + "/" + params.file,
          format: params.format,
          identity: credentials.username,
          secret: credentials.password
      };
  };

  var toHttpResource = function(name, params, credentials) {
      return {
          name: name,
          url: params.url,
          format: params.format,
          identity: credentials.username,
          secret: credentials.password
      };
  };

  var toLocalResource = function(name, params, credentials) {
    return {
      name: name,
      url: "file://" + params.path,
      format: params.format
    };
  };

  var toOpalResource = function(name, params, credentials) {
      return {
          name: name,
          url: "opal+" + params.url + "/ws/files" + params.path,
          format: params.format,
          identity: null,
          secret: credentials.token
      }
  };

  var toScpResource = function(name, params, credentials) {
      var port = params.port;
      if (!port || port<=0) {
          port = 22;
      }
      var path = params.path;
      if (!path.startsWith("/")) {
          path = "/" + path;
      }
      return {
          name: name,
          url: "scp://" + params.host + ":" + port + path,
          format: params.format,
          identity: credentials.username,
          secret: credentials.password
      }
  };

  //
  // Resource factory functions by resource form type
  //
  var toResourceFactories = {
    "gridfs-gds-file": toGridfsResource,
    "http-gds-file": toHttpResource,
    "local-gds-file": toLocalResource,
    "opal-gds-file": toOpalResource,
    "scp-gds-file": toScpResource
  };

  // Check if there is a resource factory function for the requested resource form type
  if (toResourceFactories[type]) {
      return toResourceFactories[type](name, params, credentials);
  }
  return undefined;
};
