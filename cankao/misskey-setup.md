---
title: Docker安装misskey
published: 2024-09-06
description: "Docker安装misskey"
image: ""
tags: ["Docker", "Linux", "Web项目"]
category: "指南"
draft: false
---

## 配置 docker-compose.yml 文件

```yml
services:
  web:
    image: misskey/misskey
    container_name: misskey_web
    restart: always
    links:
      - db
      - redis
    #     - mcaptcha
    #     - meilisearch
    depends_on:
      db:
        condition: service_healthy
      redis:
        condition: service_healthy
    ports:
      - "3030:3030"
    networks:
      - internal_network
      - external_network
    # env_file:
    #   - .config/docker.env
    volumes:
      - ./files:/misskey/files
      - ./.config:/misskey/.config:ro

  redis:
    restart: always
    container_name: misskey_redis
    image: redis:7-alpine
    networks:
      - internal_network
    volumes:
      - ./redis:/data
    healthcheck:
      test: "redis-cli ping"
      interval: 5s
      retries: 20

  db:
    restart: always
    image: postgres:15-alpine
    container_name: misskey_db
    networks:
      - internal_network
    env_file:
      - .config/docker.env
    volumes:
      - ./db:/var/lib/postgresql/data
    healthcheck:
      test: "pg_isready -U $$POSTGRES_USER -d $$POSTGRES_DB"
      interval: 5s
      retries: 20

#  mcaptcha:
#    restart: always
#    image: mcaptcha/mcaptcha:latest
#    networks:
#      internal_network:
#      external_network:
#        aliases:
#          - localhost
#    ports:
#      - 7493:7493
#    env_file:
#      - .config/docker.env
#    environment:
#      PORT: 7493
#      MCAPTCHA_redis_URL: "redis://mcaptcha_redis/"
#    depends_on:
#      db:
#        condition: service_healthy
#      mcaptcha_redis:
#        condition: service_healthy
#
#  mcaptcha_redis:
#    image: mcaptcha/cache:latest
#    networks:
#      - internal_network
#    healthcheck:
#      test: "redis-cli ping"
#      interval: 5s
#      retries: 20

#  meilisearch:
#    restart: always
#    image: getmeili/meilisearch:v1.3.4
#    environment:
#      - MEILI_NO_ANALYTICS=true
#      - MEILI_ENV=production
#    env_file:
#      - .config/meilisearch.env
#    networks:
#      - internal_network
#    volumes:
#      - ./meili_data:/meili_data

networks:
  internal_network:
    internal: true
  external_network:
```

## 配置文件

> 创建文件夹和文件
> .config/default.yml
> .config/docker.env

### .config/default.yml

```yml
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Misskey configuration
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

#   ┌─────┐
#───┘ URL └─────────────────────────────────────────────────────

# Final accessible URL seen by a user.
# You can set url from an environment variable instead.
url: https://example.tld/

# ONCE YOU HAVE STARTED THE INSTANCE, DO NOT CHANGE THE
# URL SETTINGS AFTER THAT!

#   ┌───────────────────────┐
#───┘ Port and TLS settings └───────────────────────────────────

#
# Misskey requires a reverse proxy to support HTTPS connections.
#
#                 +----- https://example.tld/ ------------+
#   +------+      |+-------------+      +----------------+|
#   | User | ---> || Proxy (443) | ---> | Misskey (3000) ||
#   +------+      |+-------------+      +----------------+|
#                 +---------------------------------------+
#
#   You need to set up a reverse proxy. (e.g. nginx)
#   An encrypted connection with HTTPS is highly recommended
#   because tokens may be transferred in GET requests.

# The port that your Misskey server should listen on.
port: 3000

#   ┌──────────────────────────┐
#───┘ PostgreSQL configuration └────────────────────────────────

db:
  host: db
  port: 5432

  # Database name
  # You can set db from an environment variable instead.
  db: misskey

  # Auth
  # You can set user and pass from environment variables instead.
  user: example-misskey-user
  pass: example-misskey-pass

  # Whether disable Caching queries
  #disableCache: true

  # Extra Connection options
  #extra:
  #  ssl: true

dbReplications: false

# You can configure any number of replicas here
#dbSlaves:
#  -
#    host:
#    port:
#    db:
#    user:
#    pass:
#  -
#    host:
#    port:
#    db:
#    user:
#    pass:

#   ┌─────────────────────┐
#───┘ Redis configuration └─────────────────────────────────────

redis:
  host: redis
  port: 6379
  #family: 0  # 0=Both, 4=IPv4, 6=IPv6
  #pass: example-pass
  #prefix: example-prefix
  #db: 1

#redisForPubsub:
#  host: redis
#  port: 6379
#  #family: 0  # 0=Both, 4=IPv4, 6=IPv6
#  #pass: example-pass
#  #prefix: example-prefix
#  #db: 1

#redisForJobQueue:
#  host: redis
#  port: 6379
#  #family: 0  # 0=Both, 4=IPv4, 6=IPv6
#  #pass: example-pass
#  #prefix: example-prefix
#  #db: 1

#redisForTimelines:
#  host: redis
#  port: 6379
#  #family: 0  # 0=Both, 4=IPv4, 6=IPv6
#  #pass: example-pass
#  #prefix: example-prefix
#  #db: 1

#   ┌───────────────────────────┐
#───┘ MeiliSearch configuration └─────────────────────────────

# You can set scope to local (default value) or global
# (include notes from remote).

#meilisearch:
#  host: meilisearch
#  port: 7700
#  apiKey: ''
#  ssl: true
#  index: ''
#  scope: local

#   ┌───────────────┐
#───┘ ID generation └───────────────────────────────────────────

# You can select the ID generation method.
# You don't usually need to change this setting, but you can
# change it according to your preferences.

# Available methods:
# aid ... Short, Millisecond accuracy
# aidx ... Millisecond accuracy
# meid ... Similar to ObjectID, Millisecond accuracy
# ulid ... Millisecond accuracy
# objectid ... This is left for backward compatibility

# ONCE YOU HAVE STARTED THE INSTANCE, DO NOT CHANGE THE
# ID SETTINGS AFTER THAT!

id: "aidx"

#   ┌────────────────┐
#───┘ Error tracking └──────────────────────────────────────────

# Sentry is available for error tracking.
# See the Sentry documentation for more details on options.

#sentryForBackend:
#  enableNodeProfiling: true
#  options:
#    dsn: 'https://examplePublicKey@o0.ingest.sentry.io/0'

#sentryForFrontend:
#  options:
#    dsn: 'https://examplePublicKey@o0.ingest.sentry.io/0'

#   ┌─────────────────────┐
#───┘ Other configuration └─────────────────────────────────────

# Whether disable HSTS
#disableHsts: true

# Number of worker processes
#clusterLimit: 1

# Job concurrency per worker
# deliverJobConcurrency: 128
# inboxJobConcurrency: 16

# Job rate limiter
# deliverJobPerSec: 128
# inboxJobPerSec: 32

# Job attempts
# deliverJobMaxAttempts: 12
# inboxJobMaxAttempts: 8

# IP address family used for outgoing request (ipv4, ipv6 or dual)
#outgoingAddressFamily: ipv4

# Proxy for HTTP/HTTPS
#proxy: http://127.0.0.1:3128

proxyBypassHosts:
  - api.deepl.com
  - api-free.deepl.com
  - www.recaptcha.net
  - hcaptcha.com
  - challenges.cloudflare.com

# Proxy for SMTP/SMTPS
#proxySmtp: http://127.0.0.1:3128   # use HTTP/1.1 CONNECT
#proxySmtp: socks4://127.0.0.1:1080 # use SOCKS4
#proxySmtp: socks5://127.0.0.1:1080 # use SOCKS5

# Media Proxy
#mediaProxy: https://example.com/proxy

# Proxy remote files (default: true)
proxyRemoteFiles: true

# Sign to ActivityPub GET request (default: true)
signToActivityPubGet: true
# For security reasons, uploading attachments from the intranet is prohibited,
# but exceptions can be made from the following settings. Default value is "undefined".
# Read changelog to learn more (Improvements of 12.90.0 (2021/09/04)).
#allowedPrivateNetworks: [
#  '127.0.0.1/32'
#]

# Upload or download file size limits (bytes)
#maxFileSize: 262144000
```

### .config/docker.env

```env
# misskey settings
# MISSKEY_URL=https://example.tld/

# db settings
POSTGRES_PASSWORD=example-misskey-pass
# DATABASE_PASSWORD=${POSTGRES_PASSWORD}
POSTGRES_USER=example-misskey-user
# DATABASE_USER=${POSTGRES_USER}
POSTGRES_DB=misskey
# DATABASE_DB=${POSTGRES_DB}
DATABASE_URL="postgres://${POSTGRES_USER}:${POSTGRES_PASSWORD}@db:5432/${POSTGRES_DB}"
```

## 启动

```bash
docker compose up -d
```
