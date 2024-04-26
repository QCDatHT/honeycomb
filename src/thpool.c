//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/thpool.h>

#define POOL_LOCK(pool) (pthread_mutex_lock(&((pool)->lock)))
#define POOL_UNLOCK(pool) (pthread_mutex_unlock(&((pool)->lock)))
#define POOL_NOTIFY(pool) (pthread_cond_signal(&((pool)->notify)))
#define POOL_WAIT(pool) (pthread_cond_wait(&((pool)->notify), &((pool)->lock)))

typedef struct {
   void (*function)(void *);
   void *argument;
} thpool_task_t;

struct thpool_t {
   pthread_mutex_t lock;
   pthread_cond_t notify;
   pthread_t *threads;
   thpool_task_t *queue;
   int32_t thread_count;
   int32_t queue_size;
   int32_t head;
   int32_t tail;
   int32_t count;
   int32_t shutdown;
   int32_t started;
};

static void *thpool_thread(void *thpool);

int32_t thpool_free(thpool_t *pool);

thpool_t *thpool_create(int32_t thread_count, int32_t queue_size)
{
   thpool_t *pool;
   int32_t i;

   if (thread_count <= 0 || thread_count > MAX_THREADS || queue_size <= 0 || queue_size > MAX_QUEUE) return NULL;

   if ((pool = (thpool_t *)malloc(sizeof(thpool_t))) == NULL) {
      if (pool) thpool_free(pool);

      return NULL;
   }

   pool->thread_count = 0;
   pool->queue_size = queue_size;
   pool->head = pool->tail = pool->count = 0;
   pool->shutdown = pool->started = 0;

   pool->threads = (pthread_t *)malloc(sizeof(pthread_t) * thread_count);
   pool->queue = (thpool_task_t *)malloc(sizeof(thpool_task_t) * queue_size);

   if ((pthread_mutex_init(&(pool->lock), NULL) != 0) || (pthread_cond_init(&(pool->notify), NULL) != 0) || (pool->threads == NULL) || (pool->queue == NULL)) {
      if (pool) thpool_free(pool);

      return NULL;
   }

   for (i = 0; i < thread_count; i++) {
      if (pthread_create(&(pool->threads[i]), NULL, thpool_thread, (void *)pool) != 0) {
         thpool_destroy(pool);
         return NULL;
      }
      pool->thread_count++;
      pool->started++;
   }

   return pool;
}

int32_t thpool_add(thpool_t *pool, void (*function)(void *), void *argument)
{
   if (NULL == pool || NULL == function) return 1;

   if (POOL_LOCK(pool) != 0) return 1;

   int32_t next = (pool->tail + 1) % pool->queue_size;
   bool err = false;

   do {
      if (pool->count == pool->queue_size || pool->shutdown) {
         err = true;
         break;
      }
      pool->queue[pool->tail].function = function;
      pool->queue[pool->tail].argument = argument;
      pool->tail = next;
      pool->count += 1;

      if (POOL_NOTIFY(pool) != 0) {
         err = true;
         break;
      }
   } while (0);

   if (POOL_UNLOCK(pool) != 0) err = true;

   return err;
}

int32_t thpool_destroy(thpool_t *pool)
{
   bool err = false;
   if (pool == NULL) return 1;
   if (POOL_LOCK(pool) != 0) return 1;

   do {
      if (pool->shutdown) {
         err = true;
         break;
      }

      pool->shutdown = 1;

      if ((pthread_cond_broadcast(&(pool->notify)) != 0) || (POOL_UNLOCK(pool) != 0)) {
         err = true;
         break;
      }

      for (int32_t i = 0; i < pool->thread_count; i++) {
         if (pthread_join(pool->threads[i], NULL) != 0) err = true;
      }
   } while (0);

   if (!err) thpool_free(pool);

   return err;
}

int32_t thpool_free(thpool_t *pool)
{
   if (pool == NULL || pool->started > 0) return 1;

   if (pool->threads) {
      free(pool->threads);
      free(pool->queue);

      POOL_LOCK(pool);
      pthread_mutex_destroy(&(pool->lock));
      pthread_cond_destroy(&(pool->notify));
   }
   free(pool);
   return 0;
}

static void *thpool_thread(void *thpool)
{
   thpool_t *pool = (thpool_t *)thpool;
   thpool_task_t task;

   while (1) {
      POOL_LOCK(pool);
      while ((pool->count == 0) && (!pool->shutdown))
         POOL_WAIT(pool);

      if (pool->shutdown == 1) break;

      task.function = pool->queue[pool->head].function;
      task.argument = pool->queue[pool->head].argument;
      pool->head = (pool->head + 1) % pool->queue_size;
      pool->count -= 1;

      POOL_UNLOCK(pool);
      (*(task.function))(task.argument);
   }

   pool->started--;

   POOL_UNLOCK(pool);
   pthread_exit(NULL);
   return (NULL);
}

// Copyright (C) 2024 Simone Rodini; Lorenzo Rossi
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
