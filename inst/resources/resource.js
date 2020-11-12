var dsOmics = {
  settings: {
    "title": "Omics Resources",
    "description": "Provides some resources based on Bioconductor ecosystem.",
    "web": "https://github.com/isglobal-brge/dsOmics",
    "categories": [
      {
        "name": "gds-format",
        "title": "GDS format",
        "description": "Data are in Genomic Data Structure format (GDS)."
      },
      {
        "name": "bioc",
        "title": "Bioconductor",
        "description": "Data formats and processing tools are provided by [Bioconductor](https://bioconductor.org/)."
      },
      {
        "name": "ga4gh",
        "title": "GA4GH htsget",
        "description": "Database with support of the GA4GH htsget protocol."
      },
      {
        "name": "ega",
        "title": "EGA database",
        "description": "Access to EGA with the GA4GH htsget protocol."
      }
    ],
    "types": [
      /*
      {
        "name": "ega-htsget",
        "title": "EGA htsget data access",
        "description": "Access using the GA4GH htsget API the EGA archive",
        "tags": ["ega"],
        "parameters": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "items": [
            {
              "key": "file",
              "type": "string",
              "title": "File ID",
              "description": "Name of the file to be retrieved from EGA."
            },
            {
              "key": "reference",
              "type": "string",
              "title": "Reference",
              "description": "Chromosome to retrieve (1/chr1 depending on file encoding)"
            },
            {
              "key": "start",
              "type": "integer",
              "title": "Start",
              "description": "The start position of the range on the reference."
            },
            {
              "key": "end",
              "type": "integer",
              "title": "end",
              "description": "The end position of the range on the reference."
            },
            {
              "key": "format",
              "type": "string",
              "title": "Format",
              "description": "Available file extensions that can be converted to Genomic Data Structure (GDS) format.",
              "enum": [
                {
                  "key": "EGAhtsgetVCF",
                  "title": "VCF (Variant Call Format)"
                },
                {
                  "key": "EGAhtsgetBAM",
                  "title": "BAM (Binary version of a SAM)"
                }
              ]
            }
          ],
          "required": [
            "file", "reference", "start", "end", "format"
          ]
        },
        "credentials": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "description": "Credentials to access EGA.",
          "items": [
            {
              "key": "user",
              "type": "string",
              "title": "User",
              "description": "Authentication username to EGA"
            },
            {
              "key": "pass",
              "type": "string",
              "title": "Password",
              "description": "Authentication password to EGA"
            }
          ],
          "required": [
            "user", "pass"
          ]
        }
      },
      */
      {
        "name": "ga4gh-htsget",
        "title": "GA4GH htsget database",
        "description": "Database with the [htsget protocol](http://samtools.github.io/hts-specs/htsget.html) from GA4GH implemented.",
        "tags": ["ga4gh"],
        "parameters": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "items": [
            {
              "key": "host",
              "type": "string",
              "title": "Host",
              "description": "Remote host name or IP address of the GA4GH enabled server."
            },
            {
              "key": "sample",
              "type": "string",
              "title": "Sample ID",
              "description": "."
            },
            {
              "key": "reference",
              "type": "string",
              "title": "Reference",
              "description": "."
            },
            {
              "key": "start",
              "type": "integer",
              "title": "Start",
              "description": "."
            },
            {
              "key": "end",
              "type": "integer",
              "title": "end",
              "description": "."
            },
            {
              "key": "format",
              "type": "string",
              "title": "Format",
              "description": "Available file extensions that can be converted to Genomic Data Structure (GDS) format.",
              "enum": [
                {
                  "key": "GA4GHVCF",
                  "title": "VCF (Variant Call Format)"
                },
                {
                  "key": "GA4GHBAM",
                  "title": "BAM (Binary version of a SAM)"
                }
              ]
            }
          ],
          "required": [
            "host", "sample", "reference", "start", "end"
          ]
        }
      },
      {
        "name": "gridfs-gds-file",
        "title": "GDS data file - MongoDB GridFS",
        "description": "File resource in Genomic Data Structure (GDS) format or in a format that can be converted to GDS. The file will be downloaded from the GridFS file store of a MongoDB server.",
        "tags": ["gridfs", "data-file", "gds-format", "bioc"],
        "parameters": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "items": [
            {
              "key": "host",
              "type": "string",
              "title": "Host",
              "description": "Remote host name or IP address of the MongoDB server."
            },
            {
              "key": "port",
              "type": "integer",
              "title": "Port",
              "default": 27017,
              "description": "MongoDB port number."
            },
            {
              "key": "db",
              "type": "string",
              "title": "Database",
              "description": "MongoDB database name."
            },
            {
              "key": "file",
              "type": "string",
              "title": "File",
              "description": "File name."
            },
            {
              "key": "format",
              "type": "string",
              "title": "Format",
              "description": "Genomic Data Structure (GDS) format or in a format that can be converted to GDS.",
              "enum": [
                {
                  "key": "GDS",
                  "title": "GDS"
                },
                {
                  "key": "VCF2GDS",
                  "title": "VCF (Variant Call Format)"
                }
              ]
            }
          ],
          "required": [
            "host", "port", "db", "file", "format"
          ]
        },
        "credentials": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "description": "Credentials are optional.",
          "items": [
            {
              "key": "username",
              "type": "string",
              "title": "User name",
              "description": "Valid MongoDB user name."
            },
            {
              "key": "password",
              "type": "string",
              "title": "Password",
              "format": "password",
              "description": "The user's password."
            }
          ]
        }
      },
      {
        "name": "http-gds-file",
        "title": "GDS data file - HTTP",
        "description": "File resource in Genomic Data Structure (GDS) format or in a format that can be converted to GDS. The file will be downloaded from a HTTP server.",
        "tags": ["http", "data-file", "gds-format", "bioc"],
        "parameters": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "items": [
            {
              "key": "url",
              "type": "string",
              "title": "URL",
              "description": "Address to download the file."
            },
            {
              "key": "format",
              "type": "string",
              "title": "Format",
              "description": "Genomic Data Structure (GDS) format or in a format that can be converted to GDS.",
              "enum": [
                {
                  "key": "GDS",
                  "title": "GDS"
                },
                {
                  "key": "VCF2GDS",
                  "title": "VCF (Variant Call Format)"
                }
              ]
            }
          ],
          "required": [
            "url", "format"
          ]
        },
        "credentials": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "description": "Credentials are optional. If provided, `Basic` authorization header is applied.",
          "items": [
            {
              "key": "username",
              "type": "string",
              "title": "User name",
              "description": "Valid user name."
            },
            {
              "key": "password",
              "type": "string",
              "title": "Password",
              "format": "password",
              "description": "The user's password."
            }
          ]
        }
      },
      {
        "name": "local-gds-file",
        "title": "GDS data file - local",
        "description": "File resource in Genomic Data Structure (GDS) format or in a format that can be converted to GDS. The file is located in the R server file system.",
        "tags": ["local", "data-file", "gds-format", "bioc"],
        "parameters": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "items": [
            {
              "key": "path",
              "type": "string",
              "title": "Path",
              "description": "Path to the file."
            },
            {
              "key": "format",
              "type": "string",
              "title": "Format",
              "description": "Genomic Data Structure (GDS) format or in a format that can be converted to GDS.",
              "enum": [
                {
                  "key": "GDS",
                  "title": "GDS"
                },
                {
                  "key": "VCF2GDS",
                  "title": "VCF (Variant Call Format)"
                }
              ]
            }
          ],
          "required": [
            "path", "format"
          ]
        },
        "credentials": {
          "$schema": "http://json-schema.org/schema#",
          "description": "No credentials required: the file must be accessible from the R server."
        }
      },
      {
        "name": "opal-gds-file",
        "title": "GDS data file - Opal",
        "description": "File resource in Genomic Data Structure (GDS) format or in a format that can be converted to GDS. The file will be downloaded from the file store of a Opal server.",
        "tags": ["opal", "data-file", "gds-format", "bioc"],
        "parameters": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "items": [
            {
              "key": "url",
              "type": "string",
              "title": "URL",
              "description": "Opal server base URL."
            },
            {
              "key": "path",
              "type": "string",
              "title": "Path",
              "description": "Path to the file in the Opal server."
            },
            {
              "key": "format",
              "type": "string",
              "title": "Format",
              "description": "Genomic Data Structure (GDS) format or in a format that can be converted to GDS.",
              "enum": [
                {
                  "key": "GDS",
                  "title": "GDS"
                },
                {
                  "key": "VCF2GDS",
                  "title": "VCF (Variant Call Format)"
                }
              ]
            }
          ],
          "required": [
            "url", "path", "format"
          ]
        },
        "credentials": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "description": "Credentials are required and is a [personal API access token](http://opaldoc.obiba.org/en/latest/web-user-guide/my-profile.html#personal-access-tokens).",
          "items": [
            {
              "key": "token",
              "type": "string",
              "format": "password",
              "title": "Token",
              "description": "Personal access token, granting read access to the file."
            }
          ],
          "required": [
            "token"
          ]
        }
      },
      {
        "name": "scp-gds-file",
        "title": "GDS data file - SSH",
        "description": "File resource in Genomic Data Structure (GDS) format or in a format that can be converted to GDS. The file will be downloaded from a server accessible through SSH.",
        "tags": ["ssh", "data-file", "gds-format", "bioc"],
        "parameters": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "items": [
            {
              "key": "host",
              "type": "string",
              "title": "Host",
              "description": "Remote host name or IP address that exposes SSH entry point."
            },
            {
              "key": "port",
              "type": "integer",
              "title": "Port",
              "default": 22,
              "description": "SSH port number (default is 22)."
            },
            {
              "key": "path",
              "type": "string",
              "title": "Path",
              "description": "Path to the file in the remote server."
            },
            {
              "key": "format",
              "type": "string",
              "title": "Format",
              "description": "Genomic Data Structure (GDS) format or in a format that can be converted to GDS.",
              "enum": [
                {
                  "key": "GDS",
                  "title": "GDS"
                },
                {
                  "key": "VCF2GDS",
                  "title": "VCF (Variant Call Format)"
                }
              ]
            }
          ],
          "required": [
            "host", "path", "format"
          ]
        },
        "credentials": {
          "$schema": "http://json-schema.org/schema#",
          "type": "array",
          "items": [
            {
              "key": "username",
              "type": "string",
              "title": "User name",
              "description": "Valid user name having SSH access."
            },
            {
              "key": "password",
              "type": "string",
              "title": "Password",
              "format": "password",
              "description": "The user's password."
            }
          ],
          "required": [
            "username", "password"
          ]
        }
      }
    ]
  },
  asResource: function(type, name, params, credentials) {

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
        };
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
        };
    };
    
    var toHttpGA4GH = function(name, params){
      if(params.format == "GA4GHVCF"){
        url_complete = params.host + "/variants/" + params.sample + "?format=VCF&referenceName=" + params.reference + "&start=" + params.start + "&end=" + params.end;
      }
      else{
        url_complete = params.host + "/reads/" + params.sample + "?format=BAM&referenceName=" + params.reference + "&start=" + params.start + "&end=" + params.end;
      }
      return {
            name: name,
            url: url_complete,
            format: params.format,
            secret: credentials.token
        };
    };
    
    var toHttpEGA = function(name, params, credentials){
      return{
        name: name,
      url: "https://ega.ebi.ac.uk:8052/elixir/tickets/tickets/" + params.file + "&referenceName=" + params.reference + "&start=" + params.start + "&end=" + params.end,
      identity: credentials.user,
      secret: credentials.pass,
      format: params.format
      };
    };

    //
    // Resource factory functions by resource form type
    //
    var toResourceFactories = {
      "ga4gh-htsget": toHttpGA4GH,
      "ega-htsget": toHttpEGA,
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
  }
}
